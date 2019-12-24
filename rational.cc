#include "rational.hh"
#include "makros.hh"
#include "combinatoric.hh" //imports gcd

#include "integer_math.hh"

#include <cmath>

static int __attribute__((used))  rational_res_is_save = 1;

//#define DEBUG_OUTPUT

#ifdef USE_ASM
//#define VIA128

inline void simple_save_mul(Int64& n, const Int64 fac)
{
  __asm__ volatile ("imulq %1 \n\t"
                    "jno 1f \n\t"
                    "movl $0,_ZL20rational_res_is_save \n\t"
                    "1: \n\t"
                    : "+&a"(n) : "g"(fac) : "rdx"); //load n into rax 
}
#endif

bool rational_result_is_save()
{
  assert(sizeof(Int64) == 8);
  return (rational_res_is_save != 0);
}

void reset_rational_result_is_save()
{
  rational_res_is_save = 1;
}

Rational64::Rational64()
  : num_(0), denom_(1)
{
}

Rational64::Rational64(Int64 num, Int64 denom)
  : num_(num), denom_(denom)
{
  assert(denom_ > 0);
}

Rational64::Rational64(Int64 num) : num_(num), denom_(1) {}

Rational64 Rational64::inverse() const
{
  assert(num_ != 0);

  Int64 sign = (num_ < 0) ? -1 : 1;
  return Rational64(denom_ * sign, abs(num_));
}

void Rational64::invert()
{
  assert(num_ != 0);
  if (num_ < 0) {
    num_ = - num_;
    denom_ *= -denom_;
  }

  std::swap(num_,denom_);
}

void Rational64::negate()
{
  num_ = - num_;
}

Rational64 Rational64::square() const 
{
  //we can save both gcds as this number should be normalized
#ifndef USE_ASM
  return Rational64(num_* num_, denom_ * denom_);
#else
  Rational64 result(num_, denom_);
  simple_save_mul(result.num_,result.num_);
  simple_save_mul(result.denom_,result.denom_);
  return result;
#endif
}

void Rational64::square_this()
{
#ifndef USE_ASM
  num_ *= num_;
  denom_ *= denom_
#else
  simple_save_mul(num_,num_);
  simple_save_mul(denom_,denom_);
#endif
}

void Rational64::normalize()
{
  long cur_gcd = gcd64(abs(num_),denom_);

  if (cur_gcd != 1) {
    num_ /= cur_gcd;
    denom_ /= cur_gcd;
  }
}

bool Rational64::is_normalized() const
{
  return (gcd64(abs(num_),denom_) == 1);
}

double Rational64::toDouble() const
{
  return double(num_) / double(denom_);
}

void Rational64::operator*=(Int64 fac)
{
  Int64 cur_gcd = gcd64(abs(fac),denom_);
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

Rational64 operator*(Int64 r1, const Rational64& r2)
{
  Rational64 result = r2;

  Int64 cur_gcd = gcd64(abs(r1),result.denom_);

  if (cur_gcd != 1) {
    r1 /= cur_gcd;
    result.denom_ /= cur_gcd;
  }

#ifndef USE_ASM
  result.num_ *= r1;
#else
  simple_save_mul(result.num_, r1);
#endif

  return result;
}

Rational64 operator*(const Rational64& r1, long r2)
{
  return operator*(r2,r1);
}

#ifdef USE_ASM
inline void save_add(const Rational64& r1, const Rational64& r2, Int64 denom_gcd,
                     Int64& new_num, Int64& new_denom)
{
  Int64 sdenom1 = r1.denom_;
  Int64 sdenom2 = r2.denom_;

  if (denom_gcd != 1) {
    sdenom1 /= denom_gcd;
    sdenom2 /= denom_gcd;
  }

  //assembler implementation with local labels
  __asm__ volatile ("movq %3,%%rax          \n\t"
                    "imulq %6               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%%rbx       \n\t"
                    "movq %5,%%rax          \n\t"
                    "imulq %4               \n\t"
                    "jc 1f            \n\t"
                    "addq %%rbx,%%rax       \n\t"
                    "movq %%rax,%1          \n\t" //numerator is done
                    "jo 1f            \n\t"
                    "movq %7,%%rax          \n\t"
                    "imulq %6               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%2          \n\t" //denominator is done
                    "jmp 2f            \n\t"
                    "1:              \n\t"
                    "movl $0,%0             \n\t"
                    "2:               \n\t"
                    : "+&g"(rational_res_is_save), "=&g"(new_num), "=&g"(new_denom) //0,1,2
                    : "g"(r1.num_), "g"(sdenom1), "g"(r2.num_), "g"(sdenom2), "g"(r1.denom_) //3,4,5,6,7
                    : "rax", "rbx", "rdx", "cc"
                   );
}

inline void save_add2(const Rational64& r1, const Rational64& r2, Int64 denom_gcd,
                      Int64& new_num, Int64& new_denom)
{
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

Rational64 operator+(const Rational64& r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator+" << std::endl;
#endif

  assert(r1.denom_ > 0);
  assert(r2.denom_ > 0);

  Int64 denom_gcd = gcd64(r1.denom_,r2.denom_);

  Int64 new_num;
  Int64 new_denom;

#if defined(USE_ASM) && defined(VIA128)

  Int64 prod1_high, prod1_low, prod2_high, prod2_low;
  imul(r1.num_, r2.denom_/denom_gcd, prod1_high, prod1_low);
  imul(r2.num_, r1.denom_/denom_gcd, prod2_high, prod2_low);

  bool overflow = false;

  Int64 sum_high, sum_low;
  iadd(prod1_high, prod1_low, prod2_high, prod2_low, sum_high, sum_low, overflow);
  assert(overflow == false); //this really cannot happen here!

  bool sum_fits = (sum_high == 0 || sum_high == -1);

  Int64 denom_high, denom_low;
  imul(r1.denom_/denom_gcd, r2.denom_, denom_high, denom_low);

  bool denom_fits = (denom_high == 0); //since denom is positive!

  bool mark_overflow = false;

  if (sum_fits && denom_fits) {

    long new_gcd = gcd64(abs(sum_low), denom_low);
    if (new_gcd != 1) {
      new_num = sum_low / new_gcd;
      new_denom =  denom_low / new_gcd;
    }    
  }
  else if (sum_fits) {
        
    if (denom_high >= abs(sum_low)) {
      mark_overflow = true;
    }
    else {
    
      long new_gcd = gcd_mixed_128_64(denom_high, denom_low, sum_low); 

      if (denom_high >= abs(new_gcd)) {
        mark_overflow = true;
      }    
      else {
        new_num = sum_low / new_gcd;
        new_denom = idiv(denom_high, denom_low, new_gcd);        
      }
    }
  }
  else if (denom_fits) {
    
    Int64 abs_high == sum_high;
    Int64 abs_low = sum_low;
    if (sum_high < 0)
      ineg(sum_high, sum_low, abs_high, abs_low);
    
    if (abs_high >= denom_low) {
      mark_overflow = true;
    }
    else {
    
      Int64 new_gcd = gcd_mixed_128_64(abs_high, abs_low, denom_low); 
      
      if (abs_high >= new_gcd) {
        mark_overflow = true;
      }
      else {        
        new_num = idiv(sum_high, sum_low, new_gcd);
        new_denom = denom_low / new_gcd;
      }
    }
  }
  else 
    mark_overflow = true;
  
  if (mark_overflow) {
    
    new_num = 0;
    new_denom = 1;
    
    rational_res_is_save = 0;
    assert(false);
  }

#else 
 
#ifndef USE_ASM
  new_num = r1.num_*(r2.denom_/denom_gcd)  + r2.num_*(r1.denom_/denom_gcd) ;
  new_denom = (r1.denom_/denom_gcd)*r2.denom_;
#else 
  //save_add(r1, r2, denom_gcd, new_num, new_denom);
  save_add2(r1, r2, denom_gcd, new_num, new_denom);

  assert(rational_res_is_save != 0);
#endif

  Int64 new_gcd = gcd64(abs(new_num),new_denom);
  if (new_gcd != 1) {
    new_num /= new_gcd;
    new_denom /= new_gcd;
  }
#endif

  return Rational64(new_num,new_denom);
}

void Rational64::operator+=(Rational64 r)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator+=" << std::endl;
#endif

#if defined(USE_ASM) && defined(VIA128)

  Rational64 new_r = *this + r;
  num_ = new_r.num_;
  denom_ = new_r.denom_;

#else

  assert(denom_ > 0);
  assert(r.denom_ > 0);

  Int64 denom_gcd = gcd64(denom_,r.denom_);

#ifndef USE_ASM
  Int64 new_num = num_*(r.denom_/denom_gcd)  + r.num_*(denom_/denom_gcd) ;
  Int64 new_denom = (denom_/denom_gcd)*r.denom_;
#else
  Int64 new_num;
  Int64 new_denom;

  //save_add(*this, r, denom_gcd, new_num, new_denom);
  save_add2(*this, r, denom_gcd, new_num, new_denom);
  assert(rational_res_is_save != 0);
#endif

  Int64 new_gcd = gcd64(abs(new_num),new_denom);
  if (new_gcd != 1) {
    new_num /= new_gcd;
    new_denom /= new_gcd;
  }

  num_ = new_num;
  denom_ = new_denom;
#endif
}

#ifdef USE_ASM
inline void save_sub(const Rational64& r1, const Rational64& r2, Int64 denom_gcd,
                     Int64& new_num, Int64& new_denom)
{
  Int64 sdenom1 = r1.denom_;
  Int64 sdenom2 = r2.denom_;

  if (denom_gcd != 1) {
    sdenom1 /= denom_gcd;
    sdenom2 /= denom_gcd;
  }

  //assembler implementation
  __asm__ volatile ("movq %3,%%rax          \n\t"
                    "imulq %6               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%%rbx       \n\t"
                    "movq %5,%%rax          \n\t"
                    "imulq %4               \n\t"
                    "jc 1f            \n\t"
                    "subq %%rax,%%rbx       \n\t"
                    "movq %%rbx,%1          \n\t" //numerator is done
                    "jo 1f            \n\t"
                    "movq %7,%%rax          \n\t"
                    "imulq %6               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%2          \n\t" //denominator is done
                    "jmp 2f            \n\t"
                    "1:              \n\t"
                    "movl $0,%0             \n\t"
                    "2:               \n\t"
                    : "+&g"(rational_res_is_save), "=&g"(new_num), "=&g"(new_denom) //0,1,2
                    : "g"(r1.num_), "g"(sdenom1), "g"(r2.num_), "g"(sdenom2), "g"(r1.denom_) //3,4,5,6,7
                    : "rax", "rbx", "rdx", "cc", "memory"
                   );
}

inline void save_sub2(const Rational64& r1, const Rational64& r2, Int64 denom_gcd,
                      Int64& new_num, Int64& new_denom)
{
  Int64 sdenom1 = r1.denom_;
  Int64 sdenom2 = r2.denom_;

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

Rational64 operator-(const Rational64& r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator-" << std::endl;
#endif

  assert(r1.denom_ > 0);
  assert(r2.denom_ > 0);

  Int64 denom_gcd = gcd64(r1.denom_,r2.denom_);

#ifndef USE_ASM
  Int64 new_num = r1.num_*(r2.denom_/denom_gcd)  - r2.num_*(r1.denom_/denom_gcd) ;
  Int64 new_denom = (r1.denom_/denom_gcd)*r2.denom_;
#else
  Int64 new_num;
  Int64 new_denom;

  //save_sub(r1, r2, denom_gcd, new_num, new_denom);
  save_sub2(r1, r2, denom_gcd, new_num, new_denom);

  //std::cerr << "new_num: " << new_num << std::endl;

  assert(rational_res_is_save != 0);
#endif

  const Int64 new_gcd = gcd64(abs(new_num),new_denom);
  if (new_gcd != 1) {
    new_num /= new_gcd;
    new_denom /= new_gcd;
  }

  return Rational64(new_num,new_denom);
}

void Rational64::operator-=(Rational64 r)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator-=(" << (*this) << ", " << r << ")" << std::endl;
#endif

  assert(denom_ > 0);
  assert(r.denom_ > 0);

  //std::cerr << "this: " << (*this) << ", r: " << r << std::endl;

  Int64 denom_gcd = gcd64(denom_,r.denom_);

  //std::cerr << "denom_gcd: " << denom_gcd << std::endl;

#ifndef USE_ASM
  Int64 new_num = num_*(r.denom_/denom_gcd)  - r.num_*(denom_/denom_gcd) ;
  Int64 new_denom = (denom_/denom_gcd)*r.denom_;
#else
  Int64 new_num;
  Int64 new_denom;

  //save_sub(*this, r, denom_gcd, new_num, new_denom);
  save_sub2(*this, r, denom_gcd, new_num, new_denom);

  //std::cerr << "new_num: " << new_num << std::endl;
  //std::cerr << "new denom: " << new_denom << std::endl;
  //std::cerr << "result is save: " << rational_res_is_save << std::endl;

  assert(rational_res_is_save != 0);
#endif

  const Int64 new_gcd = gcd64(abs(new_num),new_denom);
  if (new_gcd != 1) {
    new_num /= new_gcd;
    new_denom /= new_gcd;
  }

  num_ = new_num;
  denom_ = new_denom;

  assert(denom_ >= 0);
}

#ifdef USE_ASM
inline  void save_mul(Int64 n1, Int64 d1, Int64 n2, Int64 d2, Int64& num, Int64& denom)
{
  // with local labels

  __asm__ volatile ("movq %3,%%rax          \n\t"
                    "imulq %4               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%1          \n\t" //numerator done
                    "movq %5,%%rax          \n\t"
                    "imulq %6               \n\t"
                    "movq %%rax,%2          \n\t" //denominator done
                    "jmp 2f            \n\t"
                    "1:              \n\t"
                    "movl $0,%0             \n\t"
                    "2:               \n\t"
                    : "+&g"(rational_res_is_save), "=&g"(num), "=&g"(denom)
                    : "g"(n1), "g"(n2), "g"(d1), "g"(d2)
                    : "rax", "rdx", "cc"
                   );
}

inline  void save_mul2(Int64 n1, Int64 d1, Int64 n2, Int64 d2, Int64& num, Int64& denom)
{
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

Rational64 operator*(const Rational64& r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator*" << std::endl;
#endif

  if (r1.num_ == 0 || r2.num_ == 0)
    return Rational64(0,1);

  Int64 rnum1 = r1.num_;
  Int64 rdenom2 = r2.denom_;

  const Int64 gcd1 = gcd64(abs(r1.num_),r2.denom_);
  if (gcd1 != 1) {
    rnum1 /= gcd1;
    rdenom2 /= gcd1;
  }

  Int64 rnum2 = r2.num_;
  Int64 rdenom1 = r1.denom_;

  const Int64 gcd2 = gcd64(abs(r2.num_),r1.denom_);

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

  Int64 inum,idenom;
  //save_mul(rnum1,rdenom1,rnum2,rdenom2,inum,idenom);
  save_mul2(rnum1,rdenom1,rnum2,rdenom2,inum,idenom);

  assert(rational_res_is_save != 0);

  return Rational64(inum,idenom);
#endif
}

void Rational64::operator*=(Rational64 r)
{
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

  Int64 rnum1 = num_;
  Int64 rdenom2 = r.denom_;

  const Int64 gcd1 = gcd64(abs(num_),r.denom_);
  if (gcd1 != 1) {
    rnum1 /= gcd1;
    rdenom2 /= gcd1;
  }

  Int64 rnum2 = r.num_;
  Int64 rdenom1 = denom_;
  const Int64 gcd2 = gcd64(abs(r.num_),denom_);
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

  Int64 inum,idenom;
  //save_mul(rnum1,rdenom1,rnum2,rdenom2,inum,idenom);
  save_mul2(rnum1,rdenom1,rnum2,rdenom2,inum,idenom);

  assert(rational_res_is_save != 0);

  num_ = inum;
  denom_ = idenom;
#endif
}

bool operator<(const Rational64& r1, const Rational64& r2)
{
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

#ifdef VIA128
  Int64 mul1_high, mul1_low, mul2_high, mul2_low;
  imul(r1.num_, r2.denom_, mul1_high, mul1_low);
  imul(r2.num_, r1.denom_, mul2_high, mul2_low);
  
  if (mul1_high != mul2_high)
    return (mul1_high < mul2_high);
  if (neg1 && neg2)
    return (UInt64(mul1_low) > UInt64(mul2_low));
  else
    return (UInt64(mul1_low) < UInt64(mul2_low));
#else
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);

  return (ratio1 < ratio2);
#endif
}

bool operator<=(const Rational64& r1, const Rational64& r2)
{
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

#ifdef VIA128
  Int64 mul1_high, mul1_low, mul2_high, mul2_low;
  imul(r1.num_, r2.denom_, mul1_high, mul1_low);
  imul(r2.num_, r1.denom_, mul2_high, mul2_low);
  
  if (mul1_high != mul2_high)
    return (mul1_high < mul2_high);
  if (neg1 && neg2)
    return (UInt64(mul1_low) >= UInt64(mul2_low));
  else
    return (UInt64(mul1_low) <= UInt64(mul2_low));
#else
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);
#endif

  return (ratio1 <= ratio2);
}

bool operator>(const Rational64& r1, const Rational64& r2)
{
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

#ifdef VIA128
  Int64 mul1_high, mul1_low, mul2_high, mul2_low;
  imul(r1.num_, r2.denom_, mul1_high, mul1_low);
  imul(r2.num_, r1.denom_, mul2_high, mul2_low);
  
  if (mul1_high != mul2_high)
    return (mul1_high > mul2_high);
  if (neg1 && neg2)
    return (UInt64(mul1_low) < UInt64(mul2_low));
  else
    return (UInt64(mul1_low) > UInt64(mul2_low));
#else
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);
#endif

  return (ratio1 > ratio2);
}

bool operator>=(const Rational64& r1, const Rational64& r2)
{
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

#ifdef VIA128
  Int64 mul1_high, mul1_low, mul2_high, mul2_low;
  imul(r1.num_, r2.denom_, mul1_high, mul1_low);
  imul(r2.num_, r1.denom_, mul2_high, mul2_low);
  
  if (mul1_high != mul2_high)
    return (mul1_high > mul2_high);
  if (neg1 && neg2)
    return (UInt64(mul1_low) < UInt64(mul2_low));
  else
    return (UInt64(mul1_low) > UInt64(mul2_low));
#else
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);
#endif

  return (ratio1 >= ratio2);
}

std::ostream& operator<<(std::ostream& out, const Rational64& r)
{
  out << r.num_ << "/" << r.denom_;
  return out;
}

Rational64 rabs(const Rational64& r)
{
  assert(r.denom_ > 0);
  return Rational64(abs(r.num_),r.denom_);
}

bool operator==(const Rational64& r1, const Rational64& r2)
{
  if (!r2.is_normalized())
    std::cerr << "r2: " << r2 << std::endl;

  assert(r1.is_normalized());
  assert(r2.is_normalized());
  assert(r1.denom_ > 0);
  assert(r2.denom_ > 0);

  return (r1.num_==r2.num_ && r1.denom_==r2.denom_);
}

bool operator!=(const Rational64& r1, const Rational64& r2)
{
  return !operator==(r1,r2);
}

//unary minus
Rational64 operator-(const Rational64& r)
{
  return Rational64(-r.num_,r.denom_);
}

Rational64 approx_r64(double d)
{
  double df = floor(d);

  //int sign = (d < 0) ? -1 : 1;

  if (d == df)
    return Rational64(Int64(df),1);

  double absd = fabs(d);

  if (absd >= 100000000.0)
    return Rational64(Int64(floor(df+0.5)),1);
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

Rational64 approx_sqrt(const Rational64& r)
{
  if (r.num_ < 0)
    return Rational64(-1);
  else
    return Rational64(sqrt(r.num_),sqrt(r.denom_));
}
