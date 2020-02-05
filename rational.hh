#ifndef RATIONAL_HH
#define RATIONAL_HH

#include <iostream>
#include "makros.hh"

bool rational_result_is_save();

void reset_rational_result_is_save();

class Rational64 {
public:

  typedef Int64 BaseType;

  Rational64();

  Rational64(Int64 num);

  Rational64(Int64 num, Int64 denom);

  Rational64 inverse() const;

  Rational64 rabs() const;

  Rational64 negative() const;

  Rational64 square() const;

  void negate();

  void invert();

  void square_this();

  void normalize();

  bool is_normalized() const;

  inline bool is_positive() const;

  inline bool is_negative() const;

  inline bool is_zero() const;

  inline bool is_nonzero() const;

  inline bool is_one() const;

  inline bool is_integer() const;

  Int64 toInt() const;

  double toDouble() const;

  long double toLongDouble() const;

  //TODO: operators for adding and subtracting an Int64

  void operator+=(Rational64 r);

  void operator-=(Rational64 r);

  void operator*=(Int64 fac);

  void operator*=(Rational64 r);

  Int64 num_;
  Int64 denom_;
};

Rational64 operator+(const Rational64& r1, const Rational64& r2);

Rational64 operator-(const Rational64& r1, const Rational64& r2);

Rational64 operator*(const Rational64& r1, const Rational64& r2);

Rational64 operator*(Int64 r1, const Rational64& r2);

Rational64 operator*(const Rational64& r1, Int64 r2);

//unary minus
Rational64 operator-(const Rational64& r);

bool operator<(const Rational64& r1, const Rational64& r2);

bool operator<=(const Rational64& r1, const Rational64& r2);

bool operator>(const Rational64& r1, const Rational64& r2);

bool operator>=(const Rational64& r1, const Rational64& r2);

bool operator==(const Rational64& r1, const Rational64& r2);

bool operator!=(const Rational64& r1, const Rational64& r2);

std::ostream& operator<<(std::ostream& out, const Rational64& r);

//absolute of a rational number
Rational64 rabs(const Rational64& r);

Rational64 approx_sqrt(const Rational64& r);

Rational64 approx_r64(double d);

namespace Makros {

  template<>
  inline Rational64 abs(Rational64 arg)
  {
    return rabs(arg);
  }

  template<>
  inline void unified_assign(Rational64* attr_restrict dest, const Rational64* attr_restrict source, size_t size)
  {
    memcpy(dest, source, size * sizeof(Rational64));
  }
}

/******* implementation ******/

inline bool Rational64::is_negative() const
{
  assert(denom_ > 0);
  return (num_ < 0);
}

inline bool Rational64::is_positive() const
{
  assert(denom_ > 0);
  return (num_ > 0);
}

inline bool Rational64::is_one() const
{
  assert(denom_ > 0);
  return (denom_ == 1 && num_ == 1);
}

inline bool Rational64::is_zero() const
{
  return (num_ == 0);
}

inline bool Rational64::is_nonzero() const
{
  return (num_ != 0);
}

inline bool Rational64::is_integer() const
{
  //assumes normalized state
  return (denom_ == 1);
}

#endif
