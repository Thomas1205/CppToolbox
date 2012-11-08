/*** written by Thomas Schoenemann at the University of Düsseldorf, 2012 ***/

#include "rational.hh"
#include "makros.hh"
#include "combinatoric.hh" //imports gcd

#include <cmath>

static int __attribute__((used))  rational_res_is_save = 1;

//#define DEBUG_OUTPUT

bool rational_result_is_save() {
  assert(sizeof(long) == 8);
  return (rational_res_is_save != 0);
}

void reset_rational_result_is_save() {
  rational_res_is_save = 1;
}

Rational64::Rational64()
  : num_(0), denom_(1) {
}

Rational64::Rational64(long num, long denom)
  : num_(num), denom_(denom) {
  assert(denom_ > 0);
}

Rational64 Rational64::inverse() const {

  assert(num_ != 0);

  long sign = (num_ < 0) ? -1 : 1;
  return Rational64(denom_ * sign, abs(num_));
}

void Rational64::invert() {

  assert(num_ != 0);
  if (num_ < 0) {
    num_ *= -1;
    denom_ *= -1;
  }

  std::swap(num_,denom_);
}


void Rational64::negate() {
  num_ *= -1;
}

void Rational64::normalize() {

  long cur_gcd = gcd64(abs(num_),denom_);

  if (cur_gcd != 1) {
    num_ /= cur_gcd;
    denom_ /= cur_gcd;
  }
}

bool Rational64::is_normalized() const {
  return (gcd64(abs(num_),denom_) == 1);
}

bool Rational64::is_negative() const {
  assert(denom_ > 0);
  return (num_ < 0);
}

bool Rational64::is_one() const {
  assert(denom_ > 0);
  return (denom_ == 1 && num_ == 1);
}

bool Rational64::is_zero() const {
  return (num_ == 0);
}

bool Rational64::is_nonzero() const {
  return (num_ != 0);
}

double Rational64::toDouble() const {
  return double(num_) / double(denom_);
}

#ifdef USE_ASM
inline void simple_save_mul(long& n, const long fac) {

  __asm__ volatile ("imulq %1 \n\t"
                    "jno 1f \n\t"
                    "movl $0,_ZL20rational_res_is_save \n\t"
                    "1: \n\t"
                    : "+&a"(n) : "g"(fac) : "rdx");

}
#endif

void Rational64::operator*=(long fac) {

  long cur_gcd = gcd64(abs(fac),denom_);
  if (cur_gcd != 1) {
    fac /= cur_gcd;
    denom_ /= cur_gcd;
  }

#ifndef USE_ASM
  num_ *= fac;
#else
  simple_save_mul(num_, fac);
#endif
}

#ifdef USE_ASM
inline void save_add2(const Rational64& r1, const Rational64& r2, long denom_gcd,
                      long& new_num, long& new_denom) {
  
  long sdenom1 = r1.denom_;
  long sdenom2 = r2.denom_;

  if (denom_gcd != 1) {
    sdenom1 /= denom_gcd;
    sdenom2 /= denom_gcd;
  }

  //assembler implementation with local labels
  __asm__ volatile ("imulq %5               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%%rbx       \n\t"
                    "movq %4,%%rax          \n\t"
                    "imulq %3               \n\t"
                    "jc 1f            \n\t"
                    "addq %%rbx,%%rax       \n\t"
                    "movq %%rax,%0          \n\t" //numerator is done 
                    "jo 1f            \n\t"
                    "movq %6,%%rax          \n\t"
                    "imulq %5               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%1          \n\t" //denominator is done 
                    "jmp 2f            \n\t"
                    "1:              \n\t"
                    "movl $0,_ZL20rational_res_is_save        \n\t"
                    "2:               \n\t"
                    : "=&g"(new_num), "=&g"(new_denom) //0,1,
                    : "a"(r1.num_), "g"(sdenom1), "g"(r2.num_), "g"(sdenom2), "g"(r1.denom_) //2,3,4,5,6
                    : "rbx", "rdx", "cc"
                    );
}
#endif

Rational64 operator+(const Rational64& r1, const Rational64& r2) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator+" << std::endl;
#endif

  assert(r1.denom_ > 0);
  assert(r2.denom_ > 0);

  long denom_gcd = gcd64(r1.denom_,r2.denom_);

#ifndef USE_ASM
  long new_num = r1.num_*(r2.denom_/denom_gcd)  + r2.num_*(r1.denom_/denom_gcd) ;
  long new_denom = (r1.denom_/denom_gcd)*r2.denom_;
#else
  long new_num;
  long new_denom;

  save_add2(r1, r2, denom_gcd, new_num, new_denom);

  assert(rational_res_is_save != 0);
#endif

  long new_gcd = gcd64(abs(new_num),new_denom);
  if (new_gcd != 1) {
    new_num /= new_gcd;
    new_denom /= new_gcd;
  }

  return Rational64(new_num,new_denom);
}

void Rational64::operator+=(Rational64 r) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator+=" << std::endl;
#endif

  assert(denom_ > 0);
  assert(r.denom_ > 0);

  long denom_gcd = gcd64(denom_,r.denom_);

#ifndef USE_ASM
  long new_num = num_*(r.denom_/denom_gcd)  + r.num_*(denom_/denom_gcd) ;
  long new_denom = (denom_/denom_gcd)*r.denom_;
#else
  long new_num;
  long new_denom;

  save_add2(*this, r, denom_gcd, new_num, new_denom);
  assert(rational_res_is_save != 0);
#endif

  long new_gcd = gcd64(abs(new_num),new_denom);
  if (new_gcd != 1) {
    new_num /= new_gcd;
    new_denom /= new_gcd;
  }

  num_ = new_num;
  denom_ = new_denom;
}

#ifdef USE_ASM
inline void save_sub2(const Rational64& r1, const Rational64& r2, long denom_gcd,
                      long& new_num, long& new_denom) {


  long sdenom1 = r1.denom_;
  long sdenom2 = r2.denom_;

  if (denom_gcd != 1) {
    sdenom1 /= denom_gcd;
    sdenom2 /= denom_gcd;
  }

  //assembler implementation
  __asm__ volatile ("imulq %5               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%%rbx       \n\t"
                    "movq %4,%%rax          \n\t"
                    "imulq %3               \n\t"
                    "jc 1f            \n\t"
                    "subq %%rax,%%rbx       \n\t"
                    "movq %%rbx,%0          \n\t" //numerator is done 
                    "jo 1f            \n\t"
                    "movq %6,%%rax          \n\t"
                    "imulq %5               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%1          \n\t" //denominator is done 
                    "jmp 2f            \n\t"
                    "1:              \n\t"
                    "movl $0,_ZL20rational_res_is_save  \n\t"
                    "2:               \n\t"
                    : "=&g"(new_num), "=&g"(new_denom) //0,1,
                    : "a"(r1.num_), "g"(sdenom1), "g"(r2.num_), "g"(sdenom2), "g"(r1.denom_) //2,3,4,5,6
                    : "rbx", "rdx", "cc", "memory"
                    );
}
#endif

Rational64 operator-(const Rational64& r1, const Rational64& r2) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator-" << std::endl;
#endif

  assert(r1.denom_ > 0);
  assert(r2.denom_ > 0);

  long denom_gcd = gcd64(r1.denom_,r2.denom_);

#ifndef USE_ASM
  long new_num = r1.num_*(r2.denom_/denom_gcd)  - r2.num_*(r1.denom_/denom_gcd) ;
  long new_denom = (r1.denom_/denom_gcd)*r2.denom_;
#else
  long new_num;
  long new_denom;
  
  save_sub2(r1, r2, denom_gcd, new_num, new_denom);

  //std::cerr << "new_num: " << new_num << std::endl;

  assert(rational_res_is_save != 0);
#endif

  const long new_gcd = gcd64(abs(new_num),new_denom);
  if (new_gcd != 1) {
    new_num /= new_gcd;
    new_denom /= new_gcd;
  }

  return Rational64(new_num,new_denom);
}

void Rational64::operator-=(Rational64 r) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator-=(" << (*this) << ", " << r << ")" << std::endl;
#endif

  assert(denom_ > 0);
  assert(r.denom_ > 0);

  //std::cerr << "this: " << (*this) << ", r: " << r << std::endl;

  long denom_gcd = gcd64(denom_,r.denom_);

  //std::cerr << "denom_gcd: " << denom_gcd << std::endl;

#ifndef USE_ASM
  long new_num = num_*(r.denom_/denom_gcd)  - r.num_*(denom_/denom_gcd) ;
  long new_denom = (denom_/denom_gcd)*r.denom_;
#else
  long new_num;
  long new_denom;
  
  save_sub2(*this, r, denom_gcd, new_num, new_denom);


  //std::cerr << "new_num: " << new_num << std::endl;
  //std::cerr << "new denom: " << new_denom << std::endl;
  //std::cerr << "result is save: " << rational_res_is_save << std::endl;

  assert(rational_res_is_save != 0);
#endif

  const long new_gcd = gcd64(abs(new_num),new_denom);
  if (new_gcd != 1) {
    new_num /= new_gcd;
    new_denom /= new_gcd;
  }

  num_ = new_num;
  denom_ = new_denom;

  assert(denom_ >= 0);
}

#ifdef USE_ASM
inline  void save_mul2(long n1, long d1, long n2, long d2,
                      long& num, long& denom) {


  // with local labels

  __asm__ volatile ("imulq %3               \n\t"
                    "jc 1f                  \n\t"
                    "movq %%rax,%0          \n\t" //numerator done
                    "movq %4,%%rax          \n\t"
                    "imulq %5               \n\t"
                    "movq %%rax,%1          \n\t" //denominator done
                    "jmp 2f                 \n\t"
                    "1:                     \n\t"
                    "movl $0,_ZL20rational_res_is_save       \n\t"
                    "2:               \n\t"
                    : "=&g"(num), "=&g"(denom)
                    : "a"(n1), "g"(n2), "g"(d1), "g"(d2)
                    : "rdx", "cc"
                    );
}
#endif


Rational64 operator*(const Rational64& r1, const Rational64& r2) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator*" << std::endl;
#endif

  if (r1.num_ == 0 || r2.num_ == 0)
    return Rational64(0,1);

  long rnum1 = r1.num_;
  long rdenom2 = r2.denom_;

  const long gcd1 = gcd64(abs(r1.num_),r2.denom_);
  if (gcd1 != 1) {
    rnum1 /= gcd1;
    rdenom2 /= gcd1;
  }

  long rnum2 = r2.num_;
  long rdenom1 = r1.denom_;

  const long gcd2 = gcd64(abs(r2.num_),r1.denom_);

  if (gcd2 != 1) {
    rnum2 /= gcd2;
    rdenom1 /= gcd2;
  }
  // long rnum1 = r1.num_ / gcd1;
  // long rnum2 = r2.num_ / gcd2;
  // long rdenom1 = r1.denom_ / gcd2;
  // long rdenom2 = r2.denom_ / gcd1;

#ifndef USE_ASM
  return Rational64(rnum1*rnum2,rdenom1*rdenom2);
#else

  long inum,idenom;
  save_mul2(rnum1,rdenom1,rnum2,rdenom2,inum,idenom);

  assert(rational_res_is_save != 0);

  return Rational64(inum,idenom);
#endif
}


void Rational64::operator*=(Rational64 r) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator*=" << std::endl;
#endif

  if (num_ == 0) {
    assert(denom_ == 1);
    return;
  }
  if (r.num_ == 0) {
    num_ = 0;
    denom_ = 1;
    return;
  }

  long rnum1 = num_;
  long rdenom2 = r.denom_;


  const long gcd1 = gcd64(abs(num_),r.denom_);
  if (gcd1 != 1) {
    rnum1 /= gcd1;
    rdenom2 /= gcd1;
  }

  long rnum2 = r.num_;
  long rdenom1 = denom_;
  const long gcd2 = gcd64(abs(r.num_),denom_);
  if (gcd2 != 1) {
    rnum2 /= gcd2;
    rdenom1 /= gcd2;
  }

  //long rnum1 = num_ / gcd1;
  //long rnum2 = r.num_ / gcd2;
  //long rdenom1 = denom_ / gcd2;
  //long rdenom2 = r.denom_ / gcd1;

#ifndef USE_ASM
  num_ = rnum1*rnum2;
  denom_ = rdenom1*rdenom2;
#else

  long inum,idenom;
  save_mul2(rnum1,rdenom1,rnum2,rdenom2,inum,idenom);

  assert(rational_res_is_save != 0);

  num_ = inum;
  denom_ = idenom;
#endif
}

bool operator<(const Rational64& r1, const Rational64& r2) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator<" << std::endl;
#endif

  //double divisions are costly, so we try a few optimizations

  bool neg1 = r1.is_negative();
  bool neg2 = r2.is_negative();

  if (neg1 && !neg2)
    return true;
  if (neg2 && !neg1)
    return false;

  if (r1.num_ == 0)
    return (r2.num_ > 0);
 
  //this should be redundant
  //if (r2.num_ == 0)
  //  return neg1; 

  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);

  return (ratio1 < ratio2);
}

bool operator<=(const Rational64& r1, const Rational64& r2) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator<=" << std::endl;
#endif

  //double divisions are costly, so we try a few optimizations

  bool neg1 = r1.is_negative();
  bool neg2 = r2.is_negative();

  if (neg1 && !neg2)
    return true;
  if (neg2 && !neg1)
    return false;

  if (r1 == r2)
    return true;

  //this should be redundant
  //if (r2.num_ == 0)
  //  return neg1; 

  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);

  return (ratio1 <= ratio2);
}


bool operator>(const Rational64& r1, const Rational64& r2) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator>" << std::endl;
#endif

  //double divisions are costly, so we try a few optimizations

  bool neg1 = r1.is_negative();
  bool neg2 = r2.is_negative();

  if (neg1 && !neg2)
    return false;
  if (neg2 && !neg1)
    return true;

  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);

  return (ratio1 > ratio2);
}

bool operator>=(const Rational64& r1, const Rational64& r2) {

#ifdef DEBUG_OUTPUT
  std::cerr << "operator>=" << std::endl;
#endif

  //double divisions are costly, so we try a few optimizations

  bool neg1 = r1.is_negative();
  bool neg2 = r2.is_negative();

  if (neg1 && !neg2)
    return false;
  if (neg2 && !neg1)
    return true;

  if (r1 == r2)
    return true;

  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);

  return (ratio1 >= ratio2);
}


std::ostream& operator<<(std::ostream& out, const Rational64& r) {
  out << r.num_ << "/" << r.denom_;
  return out;
}


Rational64 rabs(const Rational64& r) {

  assert(r.denom_ > 0);
  return Rational64(abs(r.num_),r.denom_);
}


bool operator==(const Rational64& r1, const Rational64& r2) {

  if (!r2.is_normalized())
    std::cerr << "r2: " << r2 << std::endl;

  assert(r1.is_normalized());
  assert(r2.is_normalized());
  assert(r1.denom_ > 0);
  assert(r2.denom_ > 0);

  return (r1.num_==r2.num_ && r1.denom_==r2.denom_);
}

bool operator!=(const Rational64& r1, const Rational64& r2) {
  return !operator==(r1,r2);
}

//unary minus
Rational64 operator-(const Rational64& r) {
  return Rational64(-r.num_,r.denom_);
}


Rational64 approx_r64(double d) {

  double df = floor(d);

  if (d == df)
    return Rational64(long(df),1);

  double absd = fabs(d);

  if (absd >= 100000000.0)
    return Rational64(long(floor(df+0.5)),1);
  if (absd >= 10000.0) {

    Rational64 r(floor(10000*df+0.5),10000);
    r.normalize();
    return r;
  }
  if (absd >= 100.0) {
    Rational64 r(floor(1000000*df+0.5),1000000);
    r.normalize();
    return r;    
  }
  if (absd >= 1.0) {
    Rational64 r(floor(100000000*df+0.5),100000000);
    r.normalize();
    return r;    
  }

  Rational64 r(floor(100000000000*df+0.5),100000000000);
  r.normalize();
  return r;    
}
