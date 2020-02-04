#include "rational.hh"
#include "makros.hh"
#include "combinatoric.hh" //imports gcd

#include "integer_math.hh"

#include <cmath>

static int __attribute__((used))  rational_res_is_save = 1;

//#define DEBUG_OUTPUT
//#define CHECK_COMP

#ifdef USE_ASM
#define VIA128

/*****
**** Notes on At&T assembler syntax:
**** a/b/c/d means mapping to rax, rbx, rcx and rdx (you can specify incoming and outgoing mappings for the same register)
**** + means an operand is written and read
**** = means an operand is written
**** cc means the flag register (in the list of modified registers)
**** g means a variable can lie in memory or a register
****/

inline void simple_save_mul(Int64& n, const Int64 fac)
{  
  //note: cmov can only move memory to registers
  //note: + means that the operand is read and written, & that it is early clobber
  
  //NOTE: we only need overflow set, not rdx written -> can change to the two operand imul
  
  __asm__ volatile ("imulq %[fac] \n\t"
                    "jno 1f \n\t"
                    //"movl $0,_ZL20rational_res_is_save \n\t"
                    "movl $0, %[ratsav] \n\t"
                    "1: \n\t"
                    : "+&a"(n), [ratsav] "+m"(rational_res_is_save) : [fac] "g"(fac) : "rdx"); //load n into rax 
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

Rational64::Rational64(Int64 num) : num_(num), denom_(1) 
{
}

Rational64 Rational64::negative() const
{
  return Rational64(-num_,denom_);
}

Rational64 Rational64::rabs() const
{
  return Rational64(Makros::abs(num_),denom_);
}

Rational64 Rational64::inverse() const
{
  assert(num_ != 0);

  Int64 sign = (num_ < 0) ? -1 : 1;
  return Rational64(denom_ * sign, Makros::abs(num_));
}

void Rational64::invert()
{
  assert(num_ != 0);
  if (num_ < 0) {
    num_ = -num_;
    denom_ = -denom_;
  }

  std::swap(num_,denom_);
}

void Rational64::negate()
{
  num_ = -num_;
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
  denom_ *= denom_;
#else
  simple_save_mul(num_,num_);
  simple_save_mul(denom_,denom_);
#endif
}

void Rational64::normalize()
{
  long cur_gcd = gcd64(Makros::abs(num_),denom_);

  if (cur_gcd != 1) {
    num_ /= cur_gcd;
    denom_ /= cur_gcd;
  }
}

bool Rational64::is_normalized() const
{
  if (num_ == 0)
    return (denom_ == 1);
  return (gcd64(Makros::abs(num_),denom_) == 1 && denom_ > 0);
}

Int64 Rational64::toInt() const
{
  assert(denom_ == 1);
  return num_;
}

double Rational64::toDouble() const
{
  return double(num_) / double(denom_);
}

long double Rational64::toLongDouble() const
{
  const long double n = num_;
  const long double d = denom_;
  return n / d;
}

void Rational64::operator*=(Int64 fac)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator*=(int)" << std::endl;
#endif
  if (denom_ != 1) {
    const Int64 cur_gcd = gcd64(Makros::abs(fac),denom_);
    if (cur_gcd != 1) {
      fac /= cur_gcd;
      denom_ /= cur_gcd;
    }
  }

#ifndef USE_ASM
  num_ *= fac;
#else
  simple_save_mul(num_, fac);
#endif

  //std::cerr << "set to " << (*this) << std::endl;
  assert(denom_ > 0);
}

Rational64 operator*(Int64 r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator*(int)" << std::endl;
#endif
  Rational64 result = r2;

  if (r2.denom_ != 1) {
    const Int64 cur_gcd = gcd64(Makros::abs(r1),result.denom_);

    if (cur_gcd != 1) {
      r1 /= cur_gcd;
      result.denom_ /= cur_gcd;
    }
  }

#ifndef USE_ASM
  result.num_ *= r1;
#else
  simple_save_mul(result.num_, r1);
#endif

  //std::cerr << "result " << result << std::endl;

  assert(result.denom_ > 0);
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

  //assembler implementation with local labels, mapping r1.num_ initially to rax
  __asm__ volatile ("imulq %6               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%%rbx       \n\t"
                    "movq %5,%%rax          \n\t"
                    "imulq %4               \n\t"
                    "jc 1f            \n\t"
                    "addq %%rbx,%%rax       \n\t"
                    "movq %%rax,%0          \n\t" //numerator is done
                    "jo 1f            \n\t"
                    "jz 2f            \n\t" //skip computations if add gave zero
                    "movq %7,%%rax          \n\t"
                    "imulq %6               \n\t"
                    "jc 1f            \n\t"
                    //"movq %%rax,%1          \n\t" //denominator is done - mapping to rax in the constraints
                    "jmp 2f            \n\t"
                    "1:              \n\t"
                    //"movl $0,_ZL20rational_res_is_save        \n\t"
                    "movl $0, %2 \n\t"
                    "2:               \n\t"
                    : "=&g"(new_num), "=&a"(new_denom), "+m"(rational_res_is_save) //0,1,2
                    : "a"(r1.num_), "g"(sdenom1), "g"(r2.num_), "g"(sdenom2), "g"(r1.denom_) //3,4,5,6,7
                    : "rbx", "rdx", "cc"
                   );
}

inline void save_add_via128(const Rational64& r1, const Rational64& r2, Int64 denom_gcd,
                            Int64& new_num, Int64& new_denom)
{
  //std::cerr << "save_add_via128(" << r1 << "," << r2 << "," << denom_gcd << ")" << std::endl;

  Int64 sdenom1 = r1.denom_;
  Int64 sdenom2 = r2.denom_;

  if (denom_gcd != 1) {
    sdenom1 /= denom_gcd;
    sdenom2 /= denom_gcd;
  }

  Int64 prod1_high, prod1_low, prod2_high, prod2_low;
  imul(r1.num_, sdenom2, prod1_high, prod1_low);
  imul(r2.num_, sdenom1, prod2_high, prod2_low);

  //std::cerr << "prod1: " << prod1_high << "," << prod1_low << std::endl;
  //std::cerr << "prod2: " << prod2_high << "," << prod2_low << std::endl;

  bool overflow = false;

  Int64 sum_high, sum_low;
  iadd(prod1_high, prod1_low, prod2_high, prod2_low, sum_high, sum_low, overflow);
  //std::cerr << "sum: " << sum_high << "," << sum_low << std::endl;
  assert(overflow == false); //this really cannot happen here!

  bool sum_fits = (sum_high == 0 || sum_high == -1);

  Int64 denom_high, denom_low;
  imul(sdenom1, r2.denom_, denom_high, denom_low);

  //std::cerr << "denom: " << denom_high << "," << denom_low << std::endl;

  bool denom_fits = (denom_high == 0); //since denom is positive!

  bool mark_overflow = false;

  //std::cerr << "sum_fits: " << sum_fits << ", denom_fits: " << denom_fits << std::endl;

  if (sum_fits && denom_fits) {

    new_num = sum_low;
    if (new_num == 0) {
      new_denom = 1;
      return;
    }

    new_denom = denom_low;

    const Int64 new_gcd = gcd64(Makros::abs(sum_low), denom_low);
    if (new_gcd != 1) {
      new_num /= new_gcd;
      new_denom /= new_gcd;
    }
  }
  else if (sum_fits) {
#ifdef SAFE_MODE        
    std::cerr << "sum_fits" << std::endl;
#endif        
        
    if (denom_high >= Makros::abs(sum_low)) {
      mark_overflow = true;
    }
    else {
    
      const Int64 new_gcd = gcd_mixed_128_64(denom_high, denom_low, sum_low); 

      if (denom_high >= Makros::abs(new_gcd)) {
        mark_overflow = true;
      }    
      else {
        assert(new_gcd > 1);
        new_num = sum_low / new_gcd;
        new_denom = idiv(denom_high, denom_low, new_gcd);        
      }        
    }
  }
  else if (denom_fits) {

#ifdef SAFE_MODE
    std::cerr << "denom_fits, num: " << sum_high << "," << sum_low << " --- denom: " << denom_low << std::endl;
#endif    
    
    Int64 abs_high = sum_high;
    Int64 abs_low = sum_low;
    if (sum_high < 0) {
      ineg(sum_high, sum_low, abs_high, abs_low);
    }
    
    if (abs_high >= denom_low) {
      mark_overflow = true;
    }
    else {
    
      const Int64 new_gcd = gcd_mixed_128_64(abs_high, abs_low, denom_low); 
      //std::cerr << "gcd: " << new_gcd << std::endl;
      
      if (abs_high >= new_gcd) {
        mark_overflow = true;
      }
      else {
        assert(new_gcd > 1);
        
        new_num = idiv(sum_high, sum_low, new_gcd);
        new_denom = denom_low / new_gcd;
      }
    }
  }
  else {
#ifdef SAFE_MODE    
    std::cerr << "rien a faire" << std::endl;
#endif
    mark_overflow = true;
  }
  
  if (mark_overflow) {
    
    new_num = 0;
    new_denom = 1;

#ifdef SAFE_MODE    
    print_trace();
#endif
    
    //std::cerr << "save_add_via128(" << r1 << "," << r2 << "," << denom_gcd << ")" << std::endl;
    rational_res_is_save = 0;
    assert(false);
  }
  
  //std::cerr << "return " << new_num << "/" << new_denom << std::endl;
}
#endif

Rational64 operator+(const Rational64& r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator+" << std::endl;
#endif

  assert(r1.denom_ > 0);
  assert(r2.denom_ > 0);

  if (r1.num_ == -r2.num_ && r1.denom_ == r2.denom_)
    return Rational64(0);

  Rational64 result;

#ifdef USE_ASM
  if (r1.denom_ == r2.denom_) {
    bool overflow;
    iadd(r1.num_, r2.num_, result.num_, overflow);
    //in case of overflow use routines below
    if (!overflow) {
      result.denom_ = r1.denom_;
      if (result.denom_ != 1) {
        const Int64 new_gcd = gcd64(Makros::abs(result.num_), result.denom_);
        if (new_gcd != 1) {
          result.num_ /= new_gcd;
          result.denom_ /= new_gcd;
        }
      }        
      return result;
    }
  }
#endif

  const Int64 denom_gcd = gcd64(r1.denom_,r2.denom_);

#ifndef USE_ASM
  result.num_ = r1.num_*(r2.denom_/denom_gcd) + r2.num_*(r1.denom_/denom_gcd);
  if (result.num_ == 0) {
    result.denom_ = 1;
    return result;
  }
  result.denom_ = (r1.denom_/denom_gcd)*r2.denom_;
#else
#ifdef VIA128
  save_add_via128(r1, r2, denom_gcd, result.num_, result.denom_);  
#else  
  save_add(r1, r2, denom_gcd, result.num_, result.denom_);
#endif
  assert(rational_res_is_save != 0);
#endif

#ifndef VIA128
  if (result.num_ == 0) {
    result.denom_ = 1;
  }
  else if (result.denom_ != 1) {
    const Int64 new_gcd = gcd64(Makros::abs(result.num_),result.denom_);
    if (new_gcd != 1) {
      result.num_ /= new_gcd;
      result.denom_ /= new_gcd;
    }
  }
#endif

  assert(result.denom_ > 0);
  return result;
}

void Rational64::operator+=(Rational64 r)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator+=" << std::endl;
#endif

  assert(denom_ > 0);
  assert(r.denom_ > 0);

  if (r.is_zero())
    return;

  //std::cerr << "adding " << r << " to " << (*this) << std::endl;

  if (num_ == -r.num_ && denom_ == r.denom_) {
    num_ = 0,
    denom_ = 1;
    return;
  }

#ifdef USE_ASM
  if (denom_ == r.denom_) {
    bool overflow;
    iadd_inplace(num_, r.num_, overflow);
    //in case of overflow use routines below
    if (!overflow) {
      if (denom_ != 1) {
        const Int64 new_gcd = gcd64(Makros::abs(num_), denom_);
        if (new_gcd != 1) {
          num_ /= new_gcd;
          denom_ /= new_gcd;
        }
      }
      return;
    }
  }
#endif

  const Int64 denom_gcd = gcd64(denom_,r.denom_);

  Int64 new_num, new_denom;

#ifndef USE_ASM
  new_num = num_*(r.denom_/denom_gcd)  + r.num_*(denom_/denom_gcd);
  if (new_num == 0) {
    num_ = 0;
    denom_ = 1;
    return;
  }
  new_denom = (denom_/denom_gcd)*r.denom_;
#else

#ifdef VIA128
  save_add_via128(*this, r, denom_gcd, new_num, new_denom);
#else  
  save_add(*this, r, denom_gcd, new_num, new_denom);
#endif
  assert(rational_res_is_save != 0);
#endif

#ifndef VIA128
  if (new_num == 0) {
    num_ = 0;
    denom_ = 1;
    return;
  }

  if (new_denom != 1) {
    const Int64 new_gcd = gcd64(Makros::abs(new_num), new_denom);
    if (new_gcd != 1) {
      new_num /= new_gcd;
      new_denom /= new_gcd;
    }
  }
#endif

  num_ = new_num;
  denom_ = new_denom;

  //std::cerr << "set to " << (*this) << std::endl;
  assert(denom_ > 0);
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

  //assembler implementation, mapping r1.num_ to rax
  __asm__ volatile ("imulq %6               \n\t"
                    "jc 1f            \n\t"
                    "movq %%rax,%%rbx       \n\t"
                    "movq %5,%%rax          \n\t"
                    "imulq %4               \n\t"
                    "jc 1f            \n\t"
                    "subq %%rax,%%rbx       \n\t"
                    "movq %%rbx,%0          \n\t" //numerator is done
                    "jo 1f            \n\t"
                    "jz 2f            \n\t" //skip computations if sub gave zero
                    "movq %7,%%rax          \n\t"
                    "imulq %6               \n\t"
                    "jc 1f            \n\t"
                    //"movq %%rax,%1          \n\t" //denominator is done - mapping to rax in the constraints
                    "jmp 2f            \n\t"
                    "1:              \n\t"
                    //"movl $0,_ZL20rational_res_is_save  \n\t"
                    "movl $0, %2 \n\t"
                    "2:               \n\t"
                    : "=&g"(new_num), "=&a"(new_denom), "+m"(rational_res_is_save) //0,1,2
                    : "a"(r1.num_), "g"(sdenom1), "g"(r2.num_), "g"(sdenom2), "g"(r1.denom_) //3,4,5,6,7
                    : "rbx", "rdx", "cc"
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

  if (r1.num_ == r2.num_ && r1.denom_ == r2.denom_)
    return Rational64(0);

  Rational64 result;

#ifdef USE_ASM
  if (r1.denom_ == r2.denom_) {
    bool overflow;
    iadd(r1.num_, -r2.num_, result.num_, overflow);
    //in case of overflow use routines below
    if (!overflow && r2.num_ != MIN_LONG) { //negation fails for MIN_LONG
      result.denom_ = r1.denom_;
      if (result.denom_ != 1) {
        const Int64 new_gcd = gcd64(Makros::abs(result.num_), result.denom_);
        if (new_gcd != 1) {
          result.num_ /= new_gcd;
          result.denom_ /= new_gcd;
        }
      }
      return result;
    }
  }
#endif

  const Int64 denom_gcd = gcd64(r1.denom_,r2.denom_);

#ifndef USE_ASM
  result.num_ = r1.num_*(r2.denom_/denom_gcd) - r2.num_*(r1.denom_/denom_gcd);
  if (result.num_ == 0) {
    result.denom_ = 1;
    return result;
  }
  result.denom_ = (r1.denom_/denom_gcd)*r2.denom_;
#else

#ifdef VIA128
  assert(r2.num_ != MIN_LONG); //for MIN_LONG negation doesn't work
  save_add_via128(r1, r2.negative(), denom_gcd, result.num_, result.denom_);
#else
  save_sub(r1, r2, denom_gcd, result.num_, result.denom_);
#endif
  assert(rational_res_is_save != 0);
#endif

#ifndef VIA128
  if (result.num_ == 0) {
    result.denom_ = 1;
  }
  else if (result.denom_ != 1) {
    const Int64 new_gcd = gcd64(Makros::abs(result.num_),result.denom_);
    if (new_gcd != 1) {
      result.num_ /= new_gcd;
      result.denom_ /= new_gcd;
    }
  }
#endif

  //std::cerr << "result " << result << std::endl;
  assert(result.denom_ > 0);
  return result;
}

void Rational64::operator-=(Rational64 r)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator-=(" << (*this) << ", " << r << ")" << std::endl;
#endif

  assert(denom_ > 0);
  assert(r.denom_ > 0);

  if (r.is_zero())
    return;

  if (num_ == r.num_ && denom_ == r.denom_) {
    num_ = 0,
    denom_ = 1;
    return;
  }

#ifdef USE_ASM
  if (denom_ == r.denom_) {
    bool overflow;
    iadd_inplace(num_, -r.num_, overflow);
    //in case of overflow use routines below
    if (!overflow && r.num_ != MIN_LONG) { //negation fails for MIN_LONG
      if (denom_ != 1) {
        const Int64 new_gcd = gcd64(Makros::abs(num_), denom_);
        if (new_gcd != 1) {
          num_ /= new_gcd;
          denom_ /= new_gcd;
        }
      }
      return;
    }
  }
#endif

  //std::cerr << "this: " << (*this) << ", r: " << r << std::endl;

  const Int64 denom_gcd = gcd64(denom_,r.denom_);

  //std::cerr << "denom_gcd: " << denom_gcd << std::endl;

  Int64 new_num, new_denom;

#ifndef USE_ASM
  new_num = num_*(r.denom_/denom_gcd) - r.num_*(denom_/denom_gcd);
  if (new_num == 0) {
    num_ = 0;
    denom_ = 1;
    return;
  }
  new_denom = (denom_/denom_gcd)*r.denom_;
#else

#ifdef VIA128
  assert(r.num_ != MIN_LONG); //for MIN_LONG negation doesn't work
  save_add_via128(*this, r.negative(), denom_gcd, new_num, new_denom);
#else
  save_sub(*this, r, denom_gcd, new_num, new_denom);
  //std::cerr << "new_num: " << new_num << std::endl;
  //std::cerr << "new denom: " << new_denom << std::endl;
  //std::cerr << "result is save: " << rational_res_is_save << std::endl;
#endif
  assert(rational_res_is_save != 0);
#endif

#ifndef VIA128
  if (new_num == 0) {
    num_ = 0;
    denom_ = 1;
    return;
  }

  if (new_denom != 1) {
    const Int64 new_gcd = gcd64(Makros::abs(new_num),new_denom);
    if (new_gcd != 1) {
      new_num /= new_gcd;
      new_denom /= new_gcd;
    }
  }
#endif

  num_ = new_num;
  denom_ = new_denom;

  //std::cerr << "set to " << (*this) << std::endl;
  assert(denom_ >= 0);
}

#ifdef USE_ASM
inline void save_mul(Int64 n1, Int64 d1, Int64 n2, Int64 d2, Int64& num, Int64& denom)
{  
#if 0
  // with local labels, mapping n1 to rax
  __asm__ volatile ("imulq %4               \n\t"
                    "jo 1f                  \n\t"
                    "movq %%rax,%0          \n\t" //numerator done
                    "movq %5,%%rax          \n\t"
                    "imulq %6               \n\t"
                    "jo 1f                  \n\t"
                    "jmp 2f                 \n\t"
                    "1:                     \n\t"
                    "movl $0, %2 \n\t"
                    "2:               \n\t"
                    : "=&g"(num), "=a"(denom), "+m"(rational_res_is_save) //0,1,2
                    : "a"(n1), "g"(n2), "g"(d1), "g"(d2) //3,4,5,6
                    : "rdx", "cc"
                   );
#else
  //NOTE: we only need overflow set, not rdx written -> can change to the two operand imul
  __asm__ volatile ("imulq %%rbx, %%rax   \n\t" //dest is last
                    "jo 1f                \n\t"
                    "imulq %%rdx, %%rcx   \n\t" //dest is last
                    "jo 1f                \n\t"
                    "movl $0, %2          \n\t"
                    "2:                   \n\t"
                    : "=a"(num), "=c"(denom), "+m"(rational_res_is_save) //0,1,2
                    : "%a"(n1), "%b"(n2), "c"(d1), "d"(d2) //declare n1 and n2 commitative with the percent sign
                    : "cc");
#endif  
}
#endif

Rational64 operator*(const Rational64& r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator*(" << r1 << ", " << r2 << ")" << std::endl;
#endif

  assert(rational_res_is_save != 0);

  Rational64 result(0,1);

  if (r1.num_ == 0 || r2.num_ == 0)
    return result;

#if 1
  if (r1.denom_ == 1 && r2.denom_ == 1) {
#ifndef USE_ASM
    return Rational64(r1.num_*r2.num_,1);
#else 
    imul(r1.num_, r2.num_, result.num_, rational_res_is_save);
    return result;
#endif
  }
  if (r1.num_ == 1 && r2.num_ == 1) {
#ifndef USE_ASM
    return Rational64(r1.num_*r2.num_,1);
#else 
    imul(r1.denom_, r2.denom_, result.denom_, rational_res_is_save);
    result.num_ = 1;
    return result;
#endif
  }
#endif

  Int64 rnum1 = r1.num_;
  Int64 rdenom2 = r2.denom_;

  if (r2.denom_ != 1) {
    const Int64 gcd1 = gcd64(Makros::abs(r1.num_),r2.denom_);
    if (gcd1 != 1) {
      rnum1 /= gcd1;
      rdenom2 /= gcd1;
    }
  }

  Int64 rnum2 = r2.num_;
  Int64 rdenom1 = r1.denom_;

  if (r1.denom_ != 1) {
    const Int64 gcd2 = gcd64(Makros::abs(r2.num_),r1.denom_);
    if (gcd2 != 1) {
      rnum2 /= gcd2;
      rdenom1 /= gcd2;
    }
  }

#ifndef USE_ASM
  return Rational64(rnum1*rnum2,rdenom1*rdenom2);
#else
  save_mul(rnum1,rdenom1,rnum2,rdenom2,result.num_,result.denom_);

  assert(rational_res_is_save != 0);
  assert(result.denom_ > 0);
  //std::cerr << "return " << result << std::endl;
  return result;
#endif
}

void Rational64::operator*=(Rational64 r)
{
#ifdef DEBUG_OUTPUT
  std::cerr << "operator*=(" << r << ")" << std::endl;
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
  if (r.is_one())
    return;

#if 1
  if (denom_ == 1 && r.denom_ == 1) {
#ifndef USE_ASM
    num_ *= r.num_;
#else 
    imul_inplace(num_, r.num_, rational_res_is_save);
#endif
    return;
  }
  if (Makros::abs(num_) == 1 && r.num_ == 1) {
#ifndef USE_ASM
    denom_ *= r.denom_;
#else 
    imul_inplace(denom_, r.denom_, rational_res_is_save);
#endif
    return;
  }
#endif

  Int64 rnum1 = num_;
  Int64 rdenom2 = r.denom_;

  if (r.denom_ != 1) {
    const Int64 gcd1 = gcd64(Makros::abs(num_),r.denom_);
    if (gcd1 != 1) {
      rnum1 /= gcd1;
      rdenom2 /= gcd1;
    }
  }

  Int64 rnum2 = r.num_;
  Int64 rdenom1 = denom_;
  
  if (denom_ != 1) {
    const Int64 gcd2 = gcd64(Makros::abs(r.num_),denom_);
    if (gcd2 != 1) {
      rnum2 /= gcd2;
      rdenom1 /= gcd2;
    }
  }

#ifndef USE_ASM
  num_ = rnum1*rnum2;
  denom_ = rdenom1*rdenom2;
#else

  save_mul(rnum1,rdenom1,rnum2,rdenom2,num_,denom_);

  //std::cerr << "set to " << (*this) << std::endl;

  assert(rational_res_is_save != 0);
  assert(denom_ > 0);
#endif
}

bool operator<(const Rational64& r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  //std::cerr << "operator<" << std::endl;
#endif

  //long double divisions are costly, so we try a few optimizations

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

#ifdef CHECK_COMP  
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);
  bool approx = (ratio1 < ratio2);
#endif

  Int64 mul1_high, mul1_low, mul2_high, mul2_low;
  imul(r1.num_, r2.denom_, mul1_high, mul1_low);
  imul(r2.num_, r1.denom_, mul2_high, mul2_low);
  
  if (mul1_high != mul2_high) {
#ifdef CHECK_COMP      
    if ( (mul1_high < mul2_high) != approx ) {
      std::cerr << "differs1: < " << r1 << "," << r2 << std::endl; 
    }    
#endif    
    return (mul1_high < mul2_high);
  }
  else if (neg1 && neg2) {
    
    //-1 has the highest uint value
    
#ifdef CHECK_COMP  
    if ( (UInt64(mul1_low) < UInt64(mul2_low)) != approx ) {
      std::cerr << "differs2: < " << r1 << "," << r2 << std::endl; 
      std::cerr << "mul1: " << mul1_high << "," << mul1_low << std::endl;
      std::cerr << "mul2: " << mul2_high << "," << mul2_low << std::endl;
      
      std::cerr << "u1: " << UInt64(mul1_low) << std::endl;
      std::cerr << "u2: " << UInt64(mul2_low) << std::endl;      
    }
#endif
    return (UInt64(mul1_low) < UInt64(mul2_low));
  }
  else {
#ifdef CHECK_COMP      
    if ( (UInt64(mul1_low) < UInt64(mul2_low)) != approx ) {
      std::cerr << "differs3: < " << r1 << "," << r2 << std::endl; 
    }
#endif
    return (UInt64(mul1_low) < UInt64(mul2_low));
  }
#else
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);

  return (ratio1 < ratio2);
#endif
}

bool operator<=(const Rational64& r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  //std::cerr << "operator<=" << std::endl;
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

#ifdef CHECK_COMP  
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);
  bool approx = (ratio1 <= ratio2);
#endif

  Int64 mul1_high, mul1_low, mul2_high, mul2_low;
  imul(r1.num_, r2.denom_, mul1_high, mul1_low);
  imul(r2.num_, r1.denom_, mul2_high, mul2_low);
  
  if (mul1_high != mul2_high) {
#ifdef CHECK_COMP      
    if ( (mul1_high < mul2_high) != approx) {
       std::cerr << "differs1: <= " << r1 << ", " << r2 << std::endl;
    }        
#endif    
    return (mul1_high < mul2_high);
  }
  else if (neg1 && neg2) {
    
    //-1 has the highest uint value

#ifdef CHECK_COMP  
    if ( (UInt64(mul1_low) <= UInt64(mul2_low)) != approx) {
       std::cerr << "differs2: <= " << r1 << ", " << r2 << std::endl;
    }        
#endif
    return (UInt64(mul1_low) <= UInt64(mul2_low));
  }
  else {
#ifdef CHECK_COMP  
    if ( (UInt64(mul1_low) <= UInt64(mul2_low)) != approx) {
       std::cerr << "differs3: <= " << r1 << ", " << r2 << std::endl;
    }        
#endif
    return (UInt64(mul1_low) <= UInt64(mul2_low));
  }
#else
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);

  return (ratio1 <= ratio2);
#endif
}

bool operator>(const Rational64& r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  //std::cerr << "operator>" << std::endl;
#endif

  //double divisions are costly, so we try a few optimizations

  bool neg1 = r1.is_negative();
  bool neg2 = r2.is_negative();

  if (neg1 && !neg2)
    return false;
  if (neg2 && !neg1)
    return true;

#ifdef VIA128

#ifdef CHECK_COMP  
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);
  bool approx = (ratio1 > ratio2);
#endif

  Int64 mul1_high, mul1_low, mul2_high, mul2_low;
  imul(r1.num_, r2.denom_, mul1_high, mul1_low);
  imul(r2.num_, r1.denom_, mul2_high, mul2_low);
  
  if (mul1_high != mul2_high) {
#ifdef CHECK_COMP      
    if ( (mul1_high < mul2_high) != approx) {
       std::cerr << "differs1: > " << r1 << ", " << r2 << std::endl;
    }        
#endif    
    return (mul1_high > mul2_high);
  }
  if (neg1 && neg2) {
    
    //-1 has the highest uint value
#ifdef CHECK_COMP      
    if ( (UInt64(mul1_low) > UInt64(mul2_low)) != approx) {
       std::cerr << "differs2: > " << r1 << ", " << r2 << std::endl;
    }        
#endif    
    return (UInt64(mul1_low) > UInt64(mul2_low));
  }
  else {
#ifdef CHECK_COMP      
    if ( (UInt64(mul1_low) > UInt64(mul2_low)) != approx) {
       std::cerr << "differs3: > " << r1 << ", " << r2 << std::endl;
    }        
#endif
    return (UInt64(mul1_low) > UInt64(mul2_low));
  }
#else
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);

  return (ratio1 > ratio2);
#endif
}

bool operator>=(const Rational64& r1, const Rational64& r2)
{
#ifdef DEBUG_OUTPUT
  //std::cerr << "operator>=" << std::endl;
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

#ifdef CHECK_COMP  
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);
  bool approx = (ratio1 >= ratio2);
#endif

  Int64 mul1_high, mul1_low, mul2_high, mul2_low;
  imul(r1.num_, r2.denom_, mul1_high, mul1_low);
  imul(r2.num_, r1.denom_, mul2_high, mul2_low);
  
  if (mul1_high != mul2_high) {
#ifdef CHECK_COMP  
    if ( (mul1_high > mul2_high) != approx) {
       std::cerr << "differs1: >= " << r1 << ", " << r2 << std::endl;
    }    
#endif
    return (mul1_high > mul2_high);
  }
  if (neg1 && neg2) {
    
    //-1 has the highest uint value
#ifdef CHECK_COMP  
    if ( (UInt64(mul1_low) >= UInt64(mul2_low)) != approx) {
       std::cerr << "differs2: >= " << r1 << ", " << r2 << std::endl;
    }    
#endif
    return (UInt64(mul1_low) >= UInt64(mul2_low));
  }
  else {
#ifdef CHECK_COMP  
    if ( (UInt64(mul1_low) > UInt64(mul2_low)) != approx) {
       std::cerr << "differs3: >= " << r1 << ", " << r2 << std::endl;
    }    
#endif
    return (UInt64(mul1_low) > UInt64(mul2_low));
  }
#else
  long double ratio1 = ((long double)r1.num_) / ((long double)r1.denom_);
  long double ratio2 = ((long double)r2.num_) / ((long double)r2.denom_);

  return (ratio1 >= ratio2);
#endif
}

std::ostream& operator<<(std::ostream& out, const Rational64& r)
{
  out << r.num_ << "/" << r.denom_;
  return out;
}

Rational64 rabs(const Rational64& r)
{
  assert(r.denom_ > 0);
  return Rational64(Makros::abs(r.num_),r.denom_);
}

bool operator==(const Rational64& r1, const Rational64& r2)
{
#ifdef SAFE_MODE
  if (!r1.is_normalized()) {
    std::cerr << "r1: " << r1 << std::endl;
    std::cerr << "MAX_LONG: " << MAX_LONG << std::endl;
    std::cerr << "abs num: " << Makros::abs(r1.num_) << std::endl;
    std::cerr << "gcd: " << gcd64(Makros::abs(r1.num_),r1.denom_) << std::endl;
  }
  if (!r2.is_normalized())
    std::cerr << "r2: " << r2 << std::endl;
#endif

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
