/*** written by Thomas Schoenemann as a private person, since December 2019 ***/

#ifndef INTEGER_MATH_HH
#define INTEGER_MATH_HH

#include "makros.hh"

#ifdef USE_ASM

/****************** negate ****************/

inline void ineg(Int64 in_high, Int64 in_low, Int64& out_high, Int64& out_low) {

	out_high = ~in_high;
	out_low = ~in_low;
	
  asm volatile ("movq %[ol], %%rax \n\t"
								"incq %%rax \n\t"
								"movq %%rax, %[ol] \n\t"
								"jnc 1f \n\t"
								"movq %[oh], %%rax \n\t"
								"incq %%rax \n\t"
								"movq %%rax, %[oh] \n\t"
								"1:  \n\t"
                : [oh] "+m" (out_high), [ol] "+m" (out_low) : : "rax");
}

/****************** adds ******************/
inline void iadd(Int64 x, Int64 y, Int64& res, bool& overflow) {
	
	int o = 0;
	
  //assembler implementation with local labels
  asm volatile ("movq %[x], %%rax \n\t"
								"addq %[y], %%rax \n\t"
								"jno 1f \n\t"
								"movl $1, %[o] \n\t"
								"1: movq %%rax, %[res] \n\t"
                : [res] "+m" (res), [o] "+m" (o) : [x] "m" (x), [y] "m" (y) : "rax");
	
  overflow = (o != 0);
}

inline void iadd(Int64 x_high, Int64 x_low, Int64 y_high, Int64 y_low, Int64& res_high, Int64& res_low, bool& overflow)
{
	int o = 0;
	
  //assembler implementation with local labels
  asm volatile ("movq %[x_low], %%rax \n\t"
								"addq %[y_low], %%rax \n\t"
								"movq %%rax, %[res_low] \n\t"
								"movq %[x_high], %%rax \n\t"
								"adcq %[y_high], %%rax \n\t"
								"jno 1f \n\t"
								"movl $1, %[o] \n\t"
								"1: movq %%rax, %[res_high] \n\t"
                : [res_low] "+m" (res_low), [res_high] "+m" (res_high), [o] "+m" (o) 
								: [x_low] "m" (x_low), [x_high] "m" (x_high), [y_low] "m" (y_low), [y_high] "m" (y_high) : "rax");
	
  overflow = (o != 0);	
}

inline void uadd(UInt64 x, UInt64 y, UInt64& res, bool& overflow) {
	
	uint c = 0;
	
  //assembler implementation with local labels
  asm volatile ("movq %[x], %%rax \n\t"
								"addq %[y], %%rax \n\t"
								"jnc 1f \n\t"
								"movl $1, %[c] \n\t"
								"1: movq %%rax, %[res] \n\t"
                : [res] "+m" (res), [c] "+m" (c) : [x] "m" (x), [y] "m" (y) : "rax");
	
  overflow = (c != 0);
}


/***************** muls *******************/

inline void imul(Int64 x, Int64 y, Int64& res_high, Int64& res_low) {
	
  asm volatile ("movq %[x], %%rax \n\t"
                "imulq %[y] \n\t"
                "movq %%rdx, %[res_high] \n \t"
								"movq %%rax, %[res_low] \n\t"
                : [res_high] "+m" (res_high), [res_low] "+m" (res_low) : [x] "m" (x), [y] "m" (y) : "rdx", "rax");
	
}

inline void umul(UInt64 x, UInt64 y, UInt64& res_high, UInt64& res_low) {
	
  asm volatile ("movq %[x], %%rax \n\t"
                "mulq %[y] \n\t"
                "movq %%rdx, %[res_high] \n \t"
								"movq %%rax, %[res_low] \n\t"
                : [res_high] "+m" (res_high), [res_low] "+m" (res_low) : [x] "m" (x), [y] "m" (y) : "rdx", "rax");	
}

/***************** divs ********************/

inline Int64 idiv(Int64 x_high, Int64 x_low, Int64 y) {
	
  Int64 res;
	
  asm volatile ("movq %[x_high], %%rdx \n\t"
                "movq %[x_low], %%rax \n\t"
                "divq %[y] \n\t"
                "movq %%rax, %[res] \n\t"
                : [res] "+m" (res) : [x_high] "m" (x_high), [x_low] "m" (x_low), [y] "m" (y) : "rdx", "rax");
	
	return res;
}

#endif

#endif