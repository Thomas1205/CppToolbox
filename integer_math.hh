/*** written by Thomas Schoenemann as a private person, since December 2019 ***/

#ifndef INTEGER_MATH_HH
#define INTEGER_MATH_HH

#include "makros.hh"

#ifdef USE_ASM

/*****
**** Notes on At&T assembler syntax:
**** a/b/c/d means mapping to rax, rbx, rcx and rdx (you can specify incoming and outgoing mappings for the same register)
**** + means an operand is written and read
**** = means an operand is written
**** cc means the flag register (in the list of modified registers)
**** g means a variable can lie in memory or a register
**** In contrast to the Intel manuals, AT&T syntax is inverted. That means that the destination is the last argument.
****/

/****************** negate ****************/

inline void ineg(Int64 in_high, Int64 in_low, Int64& out_high, Int64& out_low)
{
  out_high = ~in_high;
  out_low = ~in_low;

  //assembler implementation with local labels
  asm volatile ("addq $1, %[ol] \n\t" //inc doesn't update the carry flag!
                "jnc 1f \n\t"
                "incq %[oh] \n\t"
                "1:  \n\t"
                : [oh] "+g" (out_high), [ol] "+g" (out_low) : : "cc");
}

/****************** adds ******************/

inline void iadd(Int64 x, Int64 y, Int64& res, bool& overflow)
{
  int o = 0;

  //assembler implementation with local labels
  asm volatile ("addq %[y], %%rax \n\t"
                "jno 1f \n\t"
                "movl $1, %[o] \n\t"
                "1: \n\t"
                : [res] "=a" (res), [o] "+g" (o) : [x] "%a" (x), [y] "g" (y) : "cc");
  //note: o needs +, percent declares commutative (marks also the following operand)

  overflow = (o != 0);
}

inline void iadd_inplace(Int64& x, Int64 y, bool& overflow)
{
  int o = 0;

  //assembler implementation with local labels
  asm volatile ("addq %[y], %[x] \n\t"
                "jno 1f \n\t"
                "movl $1, %[o] \n\t"
                "1: \n\t"
                : [x] "+r" (x), [o] "+g" (o) : [y] "g" (y) : "cc"); //note: o needs +

  overflow = (o != 0);
}

inline void iadd(Int64 x_high, Int64 x_low, Int64 y_high, Int64 y_low, Int64& res_high, Int64& res_low, bool& overflow)
{
  int o = 0;

  //assembler implementation with local labels
  asm volatile ("addq %[y_low], %%rax \n\t"
                "adcq %[y_high], %%rbx \n\t"
                "jno 1f \n\t"
                "movl $1, %[o] \n\t"
                "1: \n\t"
                : [res_low] "=a" (res_low), [res_high] "=b" (res_high), [o] "+g" (o) //note: o needs +
                : [x_low] "%a" (x_low), [y_low] "g" (y_low), [x_high] "b" (x_high), [y_high] "g" (y_high) : "cc");
  // percent declares commutative (marks also the following operand)

  overflow = (o != 0);
}

inline void uadd(UInt64 x, UInt64 y, UInt64& res, bool& overflow)
{
  uint c = 0;

  //assembler implementation with local labels
  asm volatile ("addq %[y], %%rax \n\t"
                "jnc 1f \n\t"
                "movl $1, %[c] \n\t"
                "1: \n\t"
                : [res] "=a" (res), [c] "+g" (c) : [x] "%a" (x), [y] "g" (y) : "cc");
  //note: c needs +, percent declares commutative (marks also the following operand)

  overflow = (c != 0);
}

/***************** muls *******************/

inline void imul(Int64 x, Int64 y, Int64& res_high, Int64& res_low)
{
  //assembler implementation
  asm volatile ("imulq %[y] \n\t"
                : [res_high] "=d" (res_high), [res_low] "=a" (res_low) : [x] "a" (x), [y] "g" (y) : "cc");

}

inline void imul(Int64 x, Int64 y, Int64& res, int& rational_res_is_save)
{
  uint o = 0;

  //NOTE: we only need overflow set, not rdx written -> use the two operand imul

  //assembler implementation with local labels
  asm volatile ("imulq %[y], %%rax \n\t"
                "jno 1f            \n\t"
                "movl $1, %[o]     \n\t"
                "1:                \n\t"
                : [res] "=a" (res), [o] "+g" (o) : [x] "%a" (x), [y] "g" (y) : "cc");
  //note: o needs +, percent declares commutative (marks also the following operand)

  if (o != 0)
    rational_res_is_save = 0;
}

inline void imul_inplace(Int64& x, Int64 y, int& rational_res_is_save)
{
  uint o = 0;

  //NOTE: we only need overflow set, not rdx written -> use the two operand imul

  //assembler implementation with local labels
  asm volatile ("imulq %[y], %[x] \n\t"
                "jno 1f            \n\t"
                "movl $1, %[o]     \n\t"
                "1:                \n\t"
                : [x] "+r" (x), [o] "+g" (o) : [y] "g" (y) : "cc");	//note: o needs +

  if (o != 0)
    rational_res_is_save = 0;
}

inline void umul(UInt64 x, UInt64 y, UInt64& res_high, UInt64& res_low)
{
  //assembler implementation
  asm volatile ("mulq %[y] \n\t"
                : [res_high] "=d" (res_high), [res_low] "=a" (res_low) : [x] "a" (x), [y] "g" (y) : "rdx", "cc");
}

/***************** divs ********************/

inline Int64 idiv(Int64 x_high, Int64 x_low, Int64 y)
{
  Int64 res;

  //assembler implementation with local labels
  asm volatile ("idivq %[y] \n\t"
                : [res] "=a" (res) : [x_high] "d" (x_high), [x_low] "a" (x_low), [y] "g" (y) : );

  return res;
}

#endif

#endif