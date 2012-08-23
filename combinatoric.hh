/*** written by Thomas Schoenemann as a private person without employment, November 2009 ***/

#ifndef COMBINATORIC_HH
#define COMBINATORIC_HH

#include "makros.hh"

uint fac(uint n);

long double ldfac(uint n);

uint choose(uint n, uint k);

long double ldchoose(uint n, uint k);

// greatest common divisor via the Euclidean algorithm
long long gcd64(unsigned long long n1, unsigned long long n2);

// greatest common divisor via the Euclidean algorithm
uint gcd(uint n1, uint n2);


#endif
