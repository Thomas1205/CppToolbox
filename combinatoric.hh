/*** first version written by Thomas Schoenemann as a private person without employment, November 2009 ***/
/*** refined at the University of Düsseldorf, Germany, 2012 ***/

#ifndef COMBINATORIC_HH
#define COMBINATORIC_HH

#include "makros.hh"
#include "vector.hh"
#include <vector>

uint fac(uint n) noexcept;

long double ldfac(uint n) noexcept;

uint choose(uint n, uint k) noexcept;

long double ldchoose(uint n, uint k) noexcept;

long double ldchoose(uint n, uint k, const Math1D::Vector<long double>& ld_fac) noexcept;

//NOTE: C++-17 provides std::gcd in <numeric>

// greatest common divisor via the Euclidean algorithm
uint gcd(uint n1, uint n2) noexcept;

// greatest common divisor via the Euclidean algorithm
inline UInt64 gcd64(UInt64 n1, UInt64 n2) noexcept;

//returns if n is a prime number (0,1,2,3 all count as prime)
template<typename T>
bool is_prime(const T n) noexcept;

//returns 1 for primes
template<typename T>
T lowest_divisor(const T n) noexcept;

template<typename T>
void get_prime_divisors(T n, std::vector<T>& divisors) noexcept;

/**** implementation ****/

inline UInt64 gcd64(UInt64 n1, UInt64 n2) noexcept
{
  if (n1 < n2)
    std::swap(n1,n2);

  while (n2 != 0) {
    const UInt64 t = n2;
    n2 = n1 % n2;
    n1 = t;
  }

  return n1;
}

inline UInt64 gcd_mixed_128_64(UInt64 n1_high, UInt64 n1_low, UInt64 n2) noexcept
{
  assert(n1_high > 0);

#ifdef USE_ASM
  if (n2 <= 1)
    return n2;

  UInt64 t = n2;
  asm volatile ("movq %[n1_high], %%rdx \n\t"
                "movq %[n1_low], %%rax \n\t"
                "divq %[n2] \n\t"
                "movq %%rdx, %[n2]"
                : [n2] "+g" (n2) : [n1_high] "g" (n1_high), [n1_low] "g" (n1_low) : "rdx", "rax");
  n1_low = t;

  while(n2 != 0) {
    t = n2;
    n2 = n1_low % n2;
    n1_low = t;
  }

  return n1_low;
#endif
}

//returns if n is a prime number (0,1,2,3 all count as prime)
template<typename T>
bool is_prime(const T n) noexcept
{
  static_assert(std::is_unsigned<T>::value, "not unsigned");
  
  if (n <= 3)
    return true;
  if ((n & 1) == 0)
    return false; //even an >= 4
  
  const T limit = sqrt(n + 0.1); //a non-prime must have a divisor <= its square root
  for (T i = 3; i <= limit; i += 2) { //even numbers are no use (only need to test primes)
    if ( (n % i) == 0)
      return false;
  }
  
  return true;
}

//returns 1 for primes
template<typename T>
T lowest_divisor(const T n) noexcept
{
  static_assert(std::is_unsigned<T>::value, "not unsigned");
  
  if (n <= 3)
    return 1;
  if ((n & 1) == 0) 
    return 2;
  
  const T limit = sqrt(n + 0.1); //a non-prime must have a divisor <= its square root
  for (T i = 3; i <= limit; i += 2) { //even numbers are no use (only need to test primes)
    if ( (n % i) == 0)
      return i;
  }
  
  return 1;  
}

template<typename T>
void get_prime_divisors(T n, std::vector<T>& divisors) noexcept
{
  static_assert(std::is_unsigned<T>::value, "not unsigned");

  divisors.clear();
  if (n <= 1)
    return;
  
  while ((n & 1) == 0) {
    divisors.push_back(2);
    n >>= 1;
  }
  
  T limit = sqrt(n + 0.1); //for numerical imprecisions
  
  T i = 3;
  while (i <= limit) {
    if ((n % i) == 0) {
      divisors.push_back(i);
      n /= i;
      limit = sqrt(n + 0.1); //for numerical imprecisions
    }
    else
      i+= 2; //even numbers are no use (only need to test primes)
  }  
  
  if (n > 1 && !divisors.empty())
    divisors.push_back(n);
}

#endif
