/*** first version written by Thomas Schoenemann as a private person without employment, November 2009 ***/
/*** refined at the University of Düsseldorf, Germany, 2012 ***/

#include "combinatoric.hh"

uint fac(uint n) {

  uint r = 1;
  for (uint k=2; k <= n; k++)
    r *= k;

  return r;
}

long double ldfac(uint n) {

  long double r = 1.0;

  for (uint k=2; k <= n; k++)
    r *= k;

  return r;  
}

uint choose(uint n, uint k) {

  assert(k <= n);

  if (k == n || k == 0)
    return 1;

  uint r = n;
  for (uint i=n-1; i > (n-k); i--)
    r *= i;

  for (uint i=2; i <= k; i++)
    r /= i;

  return r;
}

long double ldchoose(uint n, uint k) {

  assert(k <= n);

  if (k == n || k == 0)
    return 1;

  long double r = n;
  for (uint i=n-1; i > (n-k); i--)
    r *= i;

  for (uint i=2; i <= k; i++)
    r /= i;

  return r;
}

// greatest common divisor via the Euclidean algorithm
long long gcd64(unsigned long long n1, unsigned long long n2) {

  if (n1 < n2)
    std::swap(n1,n2);

  while (n2 != 0) {
    unsigned long t = n2;
    n2 = n1 % n2;
    n1 = t;
  }

  return n1;
}

// greatest common divisor via the Euclidean algorithm
uint gcd(uint n1, uint n2) {

  if (n1 < n2)
    std::swap(n1,n2);

  while (n2 != 0) {
    uint t = n2;
    n2 = n1 % n2;
    n1 = t;
  }

  return n1;
}
