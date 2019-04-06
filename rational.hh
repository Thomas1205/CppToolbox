#ifndef RATIONAL_HH
#define RATIONAL_HH

#include <iostream>

bool rational_result_is_save();

void reset_rational_result_is_save();

long gcd(unsigned long n1, unsigned long n2);

class Rational64 {
public:

  Rational64();

  Rational64(long num);

  Rational64(long num, long denom);

  Rational64 inverse() const;

  void invert();

  void normalize();

  bool is_normalized() const;
  
  bool is_negative() const;

  void negate();

  bool is_zero() const;

  bool is_nonzero() const;

  bool is_one() const;

  double toDouble() const;

  void operator+=(Rational64 r);

  void operator-=(Rational64 r);

  void operator*=(long fac);

  void operator*=(Rational64 r);
  
  long num_;
  long denom_;
};

Rational64 operator+(const Rational64& r1, const Rational64& r2);

Rational64 operator-(const Rational64& r1, const Rational64& r2);

Rational64 operator*(const Rational64& r1, const Rational64& r2);

Rational64 operator*(long r1, const Rational64& r2);

Rational64 operator*(const Rational64& r1, long r2);

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

Rational64 approx_r64(double d);

#endif
