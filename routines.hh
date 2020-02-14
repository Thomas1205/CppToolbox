/****** written by Thomas Schoenemann as a private person. moved here in February 2020 from makros.hh ****/
/****** Some basic routines, plus binary search algorithms ******/

#ifndef ROUTINES_HH
#define ROUTINES_HH

#include "makros.hh"

//Explanations for USE_SSE:
// - values 1-4 are SSE1-4 (note that x86_64 always has at least 2), 4 includes 4.1 and 4.2
// - 5 is AVX with 256 bit 
// - 6 is FMA
// - 7 is AVX2. You need it for 256 bit packed integer math

//NOTE: [Kusswurm, Modern X86 assembly language programming] recommends to use 128-bit AVX instructions instead of SSE, even though they often list one or two more arguments
//  My own experiments have shown that that does not require more space in the executable. But rigorously moving to AVX is TODO

//NOTE: [Intel® 64 and IA-32 Architectures Software Developer’s Manual Vol 2.] clarifies that use of z/y/xmm8-15 and r8-15 makes the instructions and hence the executable longer
//  TODO: stop using these regs

//  In contrast to the Intel manuals, AT&T syntax is inverted. That means that the destination is the last argument.

namespace Routines {
  
  /***************** downshift *****************/

  inline void downshift_uint_array(uint* data, const uint pos, const uint shift, const uint nData)
  {
    assert(shift <= nData);
    uint i = pos;
    const uint end = nData-shift;
#if !defined(USE_SSE) || USE_SSE < 2
    //for (; i < end; i++)
    //  data[i] = data[i+shift];
    memmove(data+pos,data+pos+shift,(end-pos+shift-1)*sizeof(uint));
#else

    //roughly the same performance as memmove

#if USE_SSE >= 5
    for (; i + 8 <= end; i += 8) 
    {
      uint* out_ptr = data+i;
      const uint* in_ptr = out_ptr + shift;

      asm __volatile__ ("vmovdqu %[inp], %%ymm7 \n\t"
                        "vmovdqu %%ymm7, %[outp] \n\t"
                        : [outp] "=m" (out_ptr[0]) : [inp] "m" (in_ptr[0]) : "ymm7", "memory");
    }
#endif

    //movdqu is SSE2
    for (; i + 4 <= end; i += 4) 
    {
      uint* out_ptr = data+i;
      const uint* in_ptr = out_ptr + shift;

      asm __volatile__ ("movdqu %[inp], %%xmm7 \n\t"
                        "movdqu %%xmm7, %[outp] \n\t"
                        : [outp] "=m" (out_ptr[0]) : [inp] "m" (in_ptr[0]) : "xmm7", "memory");
    }
    for (; i < end; i++) {
      data[i] = data[i+shift];
    }
#endif
  }

  inline void downshift_int_array(int* data, const uint pos, const uint shift, const uint nData)
  {
    assert(sizeof(int) == sizeof(uint));
    downshift_uint_array((uint*) data, pos, shift, nData);
  }

  inline void downshift_float_array(float* data, const uint pos, const uint shift, const uint nData)
  {
    assert(sizeof(int) == sizeof(uint));
    downshift_uint_array((uint*) data, pos, shift, nData);
  }

  inline void downshift_double_array(double* data, const uint pos, const uint shift, const uint nData)
  {
    uint i = pos;
    const uint end = nData-shift;
#if !defined(USE_SSE) || USE_SSE < 2
    //for (; i < end; i++)
    //  data[i] = data[i+shift];
    memmove(data+pos,data+pos+shift,(end-pos+shift-1)*sizeof(double));
#else

    //roughly the same performance as memmove

#if USE_SSE >= 5
    for (; i + 4 <= end; i += 4) {

      double* out_ptr = data+i;
      const double* in_ptr = out_ptr + shift;

      asm __volatile__ ("vmovupd %[inp], %%ymm7 \n\t"
                        "vmovupd %%ymm7, %[outp] \n\t"
                        : [outp] "=m" (out_ptr[0]) : [inp] "m" (in_ptr[0]) : "ymm7", "memory");
    }
#endif

    //movupd is SSE2
    for (; i + 2 <= end; i += 2) {

      double* out_ptr = data+i;
      const double* in_ptr = out_ptr + shift;

      asm __volatile__ ("movupd %[inp], %%xmm7 \n\t"
                        "movupd %%xmm7, %[outp] \n\t"
                        : [outp] "=m" (out_ptr[0]) : [inp] "m" (in_ptr[0]) : "xmm7", "memory");
    }
    for (; i < end; i++)
      data[i] = data[i+shift];
#endif
  }

  //routine for long double (based on memmove) is TODO

  template<typename T>
  inline void downshift_array(T* data, const uint pos, const uint shift, const uint nData)
  {
    uint i = pos;
    const uint end = nData-shift;
    for (; i < end; i++)
      data[i] = data[i+shift];
  }

  template<>
  inline void downshift_array(uint* data, const uint pos, const uint shift, const uint nData)
  {
    downshift_uint_array(data, pos, shift, nData);
  }

  template<>
  inline void downshift_array(int* data, const uint pos, const uint shift, const uint nData)
  {
    downshift_int_array(data, pos, shift, nData);
  }

  template<>
  inline void downshift_array(float* data, const uint pos, const uint shift, const uint nData)
  {
    downshift_float_array(data, pos, shift, nData);
  }

  //specialization for double is TODO

  /***************** upshift *****************/

  inline void upshift_uint_array(uint* data, const int pos, const int last, const int shift)
  {
    assert(shift > 0);
    int k = last;
#if !defined(USE_SSE) || USE_SSE < 2
    //for (; k >= pos+shift; k--)
    //  data[k] = data[k-shift];
    memmove(data+pos+shift,data+pos,(last-pos-shift+1)*sizeof(uint));
#else

    //roughly the same performance as memmove
    
#if USE_SSE >= 5
    for (; k-7 >= pos+shift; k -= 8) 
    {
      uint* out_ptr = data + k - 7;
      const uint* in_ptr = out_ptr - shift;

      asm __volatile__ ("vmovdqu %[inp], %%ymm7 \n\t"
                        "vmovdqu %%ymm7, %[outp] \n\t"
                        : [outp] "=m" (out_ptr[0]) : [inp] "m" (in_ptr[0]) : "ymm7", "memory");
    }
#endif

    //movdqu is SSE2
    for (; k-3 >= pos+shift; k -= 4) 
    {
      uint* out_ptr = data + k - 3;
      const uint* in_ptr = out_ptr - shift;

      asm __volatile__ ("movdqu %[inp], %%xmm7 \n\t"
                        "movdqu %%xmm7, %[outp] \n\t"
                        : [outp] "=m" (out_ptr[0]) : [inp] "m" (in_ptr[0]) : "xmm7", "memory");
    }
    
    for (; k >= pos+shift; k--) {
      data[k] = data[k-shift];
    }
#endif
  }

  inline void upshift_int_array(int* data, const int pos, const int last, const int shift)
  {
    assert(sizeof(int) == sizeof(uint));
    upshift_uint_array((uint*) data, pos, last, shift);
  }

  inline void upshift_float_array(float* data, const int pos, const int last, const int shift)
  {
    assert(sizeof(float) == sizeof(uint));
    upshift_uint_array((uint*) data, pos, last, shift);
  }

  inline void upshift_double_array(double* data, const int pos, const int last, const int shift)
  {
    assert(shift > 0);
    int k = last;
#if !defined(USE_SSE) || USE_SSE < 2
    //for (; k >= pos+shift; k--)
    //  data[k] = data[k-shift];    
    memmove(data+pos+shift,data+pos,(last-pos-shift+1)*sizeof(double));
#else

    //roughly the same speed as memove

#if USE_SSE >= 5
    for (; k-4-shift > pos; k -= 4) 
    {
      double* out_ptr = data + k - 3;
      const double* in_ptr = out_ptr - shift;

      asm __volatile__ ("vmovupd %[inp], %%ymm7 \n\t"
                        "vmovupd %%ymm7, %[outp] \n\t"
                        : [outp] "=m" (out_ptr[0]) : [inp] "m" (in_ptr[0]) : "ymm7", "memory");
    }
#endif

    //movupd is SSE2
    for (; k-2-shift > pos; k -= 2) 
    {
      double* out_ptr = data + k - 1;
      const double* in_ptr = out_ptr - shift;

      asm __volatile__ ("movupd %[inp], %%xmm7 \n\t"
                        "movupd %%xmm7, %[outp] \n\t"
                        : [outp] "=m" (out_ptr[0]) : [inp] "m" (in_ptr[0]) : "xmm7", "memory");
    }
    for (; k >= pos+shift; k--)
      data[k] = data[k-1];
#endif
  }

  template<typename T>
  inline void upshift_array(T* data, const int pos, const int last, const int shift)
  {
    assert(shift > 0);
    for (int k = last; k >= pos+shift; k--)
      data[k] = data[k-shift];
  }

  template<>
  inline void upshift_array(uint* data, const int pos, const int last, const int shift)
  {
    upshift_uint_array(data, pos, last, shift);
  }

  template<>
  inline void upshift_array(int* data, const int pos, const int shift, const int last)
  {
    upshift_int_array(data, pos, shift, last);
  }

  template<>
  inline void upshift_array(float* data, const int pos, const int shift, const int last)
  {
    upshift_float_array(data, pos, shift, last);
  }

  //specialization for double is TODO

  /***************** find unique *****************/

  //data should contain key at most once
  inline uint find_unique_uint(const uint* data, const uint key, const uint nData)
  {
    uint i = 0;
#if !defined(USE_SSE) || USE_SSE < 5
    for (; i < nData; i++) {
      if (data[i] == key)
        return i;
    }
#else

    //NOTE: if data is not unique, i.e. contains key more than once, this may return incorrect results
    //  if occurences are close together, it may return the sum of their positions

    //NOTE: we can go for ymm registers, if we use VEXTRACTF128 to send down the upper half

    //std::cerr << "find_unique_uint(key: " << key << ")" << std::endl;

#if USE_SSE >= 7
    //need AVX2 for 256 bit packed integers. AVX2 has vpbroadcastd 

    if (nData >= 48) {

      static const uint ind[8] = {0, 1, 2, 3, 4, 5, 6, 7, 8};

      //TODO
    }
#endif
    if (nData - i >= 12) {
      static const uint ind[4] = {0, 1, 2, 3};

      const float finc = reinterpret<const uint, const float>(4);
      const float fkey = reinterpret<const uint, const float>(key);

      asm __volatile__ ("vxorps %%xmm2, %%xmm2, %%xmm2 \n\t" // set xmm2 (the array of found positions) to zero
                        "vmovdqu %[ind], %%xmm3  \n\t" //xmm3 contains the indices
                        "vbroadcastss %[finc], %%xmm4 \n\t"
                        "vbroadcastss %[fkey], %%xmm5 \n\t"
                        : : [ind] "m" (ind[0]), [finc] "m" (finc), [fkey] "m" (fkey)
                        : "xmm2", "xmm3", "xmm4", "xmm5");

      register uint res = MAX_UINT;
      for (; i+4 <= nData; i+=4) {

        //NOTE: jumps outside of asm blocks are allowed only in asm goto, but that cannot have outputs (and local variables do not seem to get assembler names)
        // => probably best to write the loop in assembler, too

        // Assembler wishlist: horizontal min and max -> would make the uniqueness assumption superflous (not even present in AVX-512)

        asm __volatile__ ("vmovdqu %1, %%xmm0  \n\t"
                          "vpcmpeqd %%xmm5, %%xmm0, %%xmm0 \n\t" //xmm0 is overwritten with mask (all 1s on equal)
                          "pblendvb %%xmm3, %%xmm2  \n\t" //if xmm0 flags 1, the index is written
                          "vptest %%xmm0, %%xmm0 \n\t" //sets the zero flag iff xmm0 is all 0
                          "jz 1f \n\t" //jump if no equals
                          "phaddd %%xmm0, %%xmm2 \n\t" //(xmm0 is irrelevant)
                          "phaddd %%xmm0, %%xmm2 \n\t" //(xmm0 is irrelevant)
                          "vpextrd $0, %%xmm2, %0 \n\t"
                          "1: paddd %%xmm4, %%xmm3 \n\t"
                          : "+g" (res) : "m" (data[i]) : "xmm0", "xmm2", "xmm3");

        if (res != MAX_UINT)
          return res;
      }
    }

    for (; i < nData; i++) {
      if (data[i] == key)
        return i;
    }
#endif
    return MAX_UINT;
  }

  //data should contain key at most once
  inline uint find_unique_int(const int* data, const int key, const uint nData)
  {
    return find_unique_uint((const uint*) data, (const uint) key, nData);
  }

  inline uint find_unique_float(const float* data, const float key, const uint nData)
  {
    return find_unique_uint((const uint*) data, reinterpret<const uint, const float>(key), nData);
  }

  /***************** find first *****************/

  inline uint find_first_uint(const uint* data, const uint key, const uint nData)
  {
    uint i = 0;
#if !defined(USE_SSE) || USE_SSE < 5
    for (; i < nData; i++) {
      if (data[i] == key)
        return i;
    }
#else

    //NOTE: if data is not unique, i.e. contains key more than once, this may return incorrect results
    //  if occurences are close together, it may return the sum of their positions

    //NOTE: we can go for ymm registers, if we use VEXTRACTF128 to send down the upper half

    std::cerr << "find_first_uint(key: " << key << ")" << std::endl;

    if (nData >= 12) {
      static const uint ind[4] = {0, 1, 2, 3};

      const float finc = reinterpret<const uint, const float>(4); //*reinterpret_cast<const float*>(&iinc);
      const float fkey = reinterpret<const uint, const float>(key); //*reinterpret_cast<const float*>(&key);
      const float fmax = reinterpret<const uint, const float>(MAX_UINT);

      asm __volatile__ ("vbroadcastss %[fmax], %%xmm2 \n\t" // set xmm2 (the array of found positions) to MAX_UINT
                        "vmovdqu %[ind], %%xmm3  \n\t" //xmm3 contains the indices
                        "vbroadcastss %[finc], %%xmm4 \n\t"
                        "vbroadcastss %[fkey], %%xmm5 \n\t"
                        : : [ind] "m" (ind[0]), [finc] "m" (finc), [fkey] "m" (fkey), [fmax] "m" (fmax)
                        : "xmm2", "xmm3", "xmm4", "xmm5");

      register uint res = MAX_UINT;
      for (; i+4 <= nData; i+=4) {

        //NOTE: jumps outside of asm blocks are allowed only in asm goto, but that cannot have outputs (and local variables do not seem to get assembler names)
        // => probably best to write the loop in assembler, too

        // Assembler wishlist: horizontal min and max -> would save the lengthy manual code

        asm __volatile__ ("vmovdqu %[dat], %%xmm0  \n\t"
                          "vpcmpeqd %%xmm5, %%xmm0, %%xmm0 \n\t" //xmm0 is overwritten with mask (all 1s on equal)
                          "pblendvb %%xmm3, %%xmm2  \n\t" //if xmm0 flags 1, the index is written
                          "vptest %%xmm0, %%xmm0 \n\t" //sets the zero flag iff xmm0 is all 0
                          "jz 1f \n\t" //jump if no equals
                          //manual phmin
                          "vmovhlps %%xmm2, %%xmm1, %%xmm1 \n\t" //move high two ints of xmm2 to low in xmm1
                          //"vpsrldq $8, %%xmm2, %%xmm1 \n\t" //move high two ints of xmm2 to low in xmm1: dq is byte shift
                          "pminud %%xmm1, %%xmm2   \n\t"
                          "vpsrldq $4, %%xmm2, %%xmm1 \n\t"
                          "pminud %%xmm1, %%xmm2   \n\t"
                          "vpextrd $0, %%xmm2, %0 \n\t"
                          "1: paddd %%xmm4, %%xmm3 \n\t"
                          : "+g" (res) : [dat] "m" (data[i]) : "xmm0", "xmm1", "xmm2", "xmm3");

        if (res != MAX_UINT)
          return res;
      }
    }

    for (; i < nData; i++) {
      if (data[i] == key)
        return i;
    }
#endif
    return MAX_UINT;
  }

  /******** min, max, min+arg_min, max+arg_max *******/

  inline float max(const float_A16* data, size_t nData)
  {
    float max_val=MIN_FLOAT;
    float cur_datum;
    size_t i = 0;

//#if !defined(USE_SSE) || USE_SSE < 2
#if 1 // g++ 4.8.5 uses avx instructions automatically
    for (; i < nData; i++) {
      cur_datum = data[i];
      max_val = std::max(max_val,cur_datum);
    }
#else
    //movaps is part of SSE2

    float tmp[4] = {MIN_FLOAT,MIN_FLOAT,MIN_FLOAT,MIN_FLOAT}; //reused as output, static not useful

    //maxps can take an unaligned mem arg!
    asm __volatile__ ("movaps %[tmp], %%xmm6" : : [tmp] "m" (tmp[0]) : "xmm6");
    for (; (i+4) <= nData; i += 4) {
      asm __volatile__ ("movaps %[fptr], %%xmm7\n\t"
                        "maxps %%xmm7, %%xmm6" : : [fptr] "m" (data[i]) : "xmm6", "xmm7");

    }
    asm __volatile__ ("movups %%xmm6, %[tmp]" : [tmp] "=m" (tmp[0]) : : "memory");
    for (k=0; k < 4; k++)
      max_val = std::max(max_val,tmp[k]);

    for (; i < nData; i++) {
      cur_datum = data[i];
      if (cur_datum > max_val)
        max_val = cur_datum;
    }
#endif

    return max_val;
  }

  inline float min(const float_A16* data, size_t nData)
  {
    float min_val=MAX_FLOAT;
    float cur_datum;
    size_t i = 0;

    //#if !defined(USE_SSE) || USE_SSE < 2
#if 1 // g++ 4.8.5 uses avx instructions automatically
    for (; i < nData; i++) {
      cur_datum = data[i];
      min_val = std::min(min_val,cur_datum);
    }
#else
    //movaps is part of SSE2

    float tmp[4] = {MAX_FLOAT,MAX_FLOAT,MAX_FLOAT,MAX_FLOAT}; //reused as output, static not useful

    //minps can take an unaligned mem arg!
    asm __volatile__ ("movaps %[tmp], %%xmm6" : : [tmp] "m" (tmp[0]) : "xmm6");
    for (; (i+4) <= nData; i += 4) {
      asm __volatile__ ("movaps %[fptr], %%xmm7 \n\t"
                        "minps %%xmm7, %%xmm6 \n\t" : : [fptr] "m" (data[i]) : "xmm6", "xmm7");
    }
    asm __volatile__ ("movups %%xmm6, %[tmp]" : [tmp] "=m" (tmp[0]) :  : "memory");
    for (k=0; k < 4; k++)
      min_val = std::min(min_val,tmp[k]);

    for (; i < nData; i++) {
      cur_datum = data[i];
      if (cur_datum < min_val)
        min_val = cur_datum;
    }
#endif

    return min_val;
  }

  inline void find_max_and_argmax(const float_A16* data, const size_t nData, float& max_val, size_t& arg_max)
  {
    max_val = MIN_FLOAT;
    arg_max = MAX_UINT;

    assertAligned16(data);

#if !defined(USE_SSE) || USE_SSE < 4

    if (nData > 0) {
      const float* ptr = std::max_element(data,data+nData);
      max_val = *ptr;
      arg_max = ptr - data;
    }

    // for (i=0; i < nData; i++) {
    //   float cur_val = data[i];

    //   if (cur_val > max_val) {
    //     max_val = cur_val;
    //     arg_max = i;
    //   }
    // }
#else    

    size_t i = 0;
    float cur_val;
    
#if USE_SSE >= 5

    //use AVX - align16 is no good, need align32 for aligned moves

    if (nData >= 16) {
      const float val = MIN_FLOAT;
      //const uint one = 1;
      const float inc = reinterpret<const uint, const float>(1); //*reinterpret_cast<const float*>(&one);
      
      //TODO: check out VPBROADCASTD (but it's AVX2)

      asm __volatile__ ("vbroadcastss %[tmp], %%ymm6 \n\t" //ymm6 is max register
                        "vxorps %%ymm5, %%ymm5, %%ymm5 \n\t" //sets ymm5 (= argmax) to zero
                        "vbroadcastss %[itemp], %%ymm4 \n\t" //ymm4 is increment register
                        "vxorps %%ymm3, %%ymm3, %%ymm3 \n\t" //sets ymm3 (= current set index) to zero
                        : : [tmp] "m" (val), [itemp] "m" (inc) : "ymm3", "ymm4", "ymm5", "ymm6");

      for (; (i+8) <= nData; i += 8) {

        asm __volatile__ ("vmovups %[fptr], %%ymm7 \n\t"
                          "vcmpnleps %%ymm6, %%ymm7, %%ymm0 \n\t"
                          "vblendvps %%ymm0, %%ymm7, %%ymm6, %%ymm6 \n\t" //destination is last
                          "vblendvps %%ymm0, %%ymm3, %%ymm5, %%ymm5 \n\t" //destination is last
                          "vpaddd %%ymm4, %%ymm3, %%ymm3 \n\t" //destination is last
                          : : [fptr] "m" (data[i]) : "ymm0", "ymm3", "ymm5", "ymm6", "ymm7");
      }

      float tmp[8];
      uint itemp[8];

      asm __volatile__ ("vmovups %%ymm6, %[tmp] \n\t"
                        "vmovups %%ymm5, %[itemp]"
                        : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "memory");

      for (uint k=0; k < 8; k++) {
        cur_val = tmp[k];
        if (cur_val > max_val) {
          max_val = cur_val;
          arg_max = (itemp[k] << 3) + k; //8*itemp[k] + k;
        }
      }
      assert(i == nData - (nData % 8));
    }

#else
    //blendvps is part of SSE4

    if (nData >= 8) {
      assert(nData <= 17179869183);

      float tmp[4] = {MIN_FLOAT,MIN_FLOAT,MIN_FLOAT,MIN_FLOAT}; //reused as output, static not useful
      uint itemp[4] = {1,1,1,1}; //increment array, reused as output, static not useful
      assert(sizeof(uint) == 4);

      asm __volatile__ ("movups %[tmp], %%xmm6 \n\t" //xmm6 is max register
                        "xorps %%xmm5, %%xmm5 \n\t" //sets xmm5 (= argmax) to zero
                        "movups %[itemp], %%xmm4 \n\t" //xmm4 is increment register
                        "xorps %%xmm3, %%xmm3 \n\t" //sets xmm3 (= current set index) to zero
                        : : [tmp] "m" (tmp[0]), [itemp] "m" (itemp[0]) : "xmm3", "xmm4", "xmm5", "xmm6");

      for (; (i+4) <= nData; i += 4) {

        asm __volatile__ ("movaps %[fptr], %%xmm7 \n\t"
                          "movaps %%xmm7, %%xmm0 \n\t"
                          "cmpnleps %%xmm6, %%xmm0 \n\t"
                          "blendvps %%xmm7, %%xmm6 \n\t"
                          "blendvps %%xmm3, %%xmm5 \n\t"
                          "paddd %%xmm4, %%xmm3 \n\t"
                          : : [fptr] "m" (data[i]) : "xmm0", "xmm3", "xmm5", "xmm6", "xmm7");
      }

      asm __volatile__ ("movups %%xmm6, %[tmp] \n\t"
                        "movups %%xmm5, %[itemp]"
                        : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "memory");

      for (uint k=0; k < 4; k++) {
        cur_val = tmp[k];
        if (cur_val > max_val) {
          max_val = cur_val;
          arg_max = 4*itemp[k] + k;
        }
      }
      assert(i == nData - (nData % 4));
    }
#endif

    for (; i < nData; i++) {
      cur_val = data[i];
      if (cur_val > max_val) {
        max_val = cur_val;
        arg_max = i;
      }
    }
#endif
  }

  inline void find_max_and_argmax(const double_A16* data, const size_t nData, double& max_val, size_t& arg_max)
  {
    max_val = MIN_DOUBLE;
    arg_max = MAX_UINT;

    assertAligned16(data);

#if !defined(USE_SSE) || USE_SSE < 4
//#if 1
    if (nData > 0) {
      const double* ptr = std::max_element(data,data+nData);
      max_val = *ptr;
      arg_max = ptr - data;
    }

    // for (i=0; i < nData; i++) {
    //   double cur_val = data[i];

    //   if (cur_val > max_val) {
    //     max_val = cur_val;
    //     arg_max = i;
    //   }
    // }
#else    

    size_t i = 0;
    double cur_val;
    
#if USE_SSE >= 5

    //use AVX  - align16 is no good, need align32 for aligned moves

    if (nData >= 8) {
      assert(sizeof(size_t) == 8);

      const double val = MIN_DOUBLE;
      const size_t one = 1;
      const double inc = reinterpret<const size_t, const double>(1); //*reinterpret_cast<const double*>(&one);

      asm __volatile__ ("vbroadcastsd %[tmp], %%ymm6 \n\t" //ymm6 is max register
                        "vxorpd %%ymm5, %%ymm5, %%ymm5 \n\t" //sets ymm5 (= argmax) to zero
                        "vbroadcastsd %[itemp], %%ymm4 \n\t" //ymm4 is increment register
                        "vxorpd %%ymm3, %%ymm3, %%ymm3 \n\t" //sets ymm3 (= current set index) to zero
                        : : [tmp] "m" (val), [itemp] "m" (inc) : "ymm3", "ymm4", "ymm5", "ymm6");

      for (; (i+4) <= nData; i += 4) {

        asm __volatile__ ("vmovupd %[dptr], %%ymm7 \n\t"
                          "vcmpnlepd %%ymm6, %%ymm7, %%ymm0 \n\t"
                          "vblendvpd %%ymm0, %%ymm7, %%ymm6, %%ymm6 \n\t" //destination is last
                          "vblendvpd %%ymm0, %%ymm3, %%ymm5, %%ymm5 \n\t" //destination is last
                          "vpaddd %%ymm4, %%ymm3, %%ymm3 \n\t" //destination is last
                          : : [dptr] "m" (data[i]) : "ymm0", "ymm3", "ymm5", "ymm6", "ymm7");
      }

      double tmp[4];
      size_t itemp[4];

      asm __volatile__ ("vmovups %%ymm6, %[tmp] \n\t"
                        "vmovups %%ymm5, %[itemp]"
                        : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "memory");

      for (uint k=0; k < 4; k++) {
        cur_val = tmp[k];
        //std::cerr << "cur val: " << cur_val << std::endl;
        if (cur_val > max_val) {
          max_val = cur_val;
          arg_max = (itemp[k] << 2) + k; //4*itemp[k] + k;
        }
      }
      //std::cerr << "minval: " << min_val << std::endl;
      assert(i == nData - (nData % 4));
    }

#else

    if (nData >= 4) {
      assert(nData < 8589934592);

      //note: broadcast is better implemented by movddup (SSE3). But this code is no longer improved
      double tmp[2] = {MIN_DOUBLE,MIN_DOUBLE}; //reused as output, static not useful
      uint itemp[4] = {0,1,0,1}; //reused as output, static not useful
      assert(sizeof(uint) == 4);

      asm __volatile__ ("movupd %[tmp], %%xmm6 \n\t"
                        "xorpd %%xmm5, %%xmm5 \n\t" //sets xmm5 (= argmax) to zero
                        "movupd %[itemp], %%xmm4 \n\t"
                        "xorpd %%xmm3, %%xmm3 \n\t" //sets xmm3 (= current set index)
                        : : [tmp] "m" (tmp[0]), [itemp] "m" (itemp[0]) :  "xmm3", "xmm4", "xmm5", "xmm6");

      for (; (i+2) <= nData; i += 2) {

        asm __volatile__ ("movapd %[dptr], %%xmm7 \n\t"
                          "movapd %%xmm7, %%xmm0 \n\t"
                          "cmpnlepd %%xmm6, %%xmm0 \n\t"
                          "blendvpd %%xmm7, %%xmm6 \n\t"
                          "blendvpd %%xmm3, %%xmm5 \n\t"
                          "paddd %%xmm4, %%xmm3 \n\t"
                          : : [dptr] "m" (data[i]) : "xmm0", "xmm3", "xmm5", "xmm6", "xmm7");
      }

      asm __volatile__ ("movupd %%xmm6, %[tmp] \n\t"
                        "movupd %%xmm5, %[itemp] \n\t"
                        : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "memory");

      assert(itemp[0] == 0);
      assert(itemp[2] == 0);

      for (uint k=0; k < 2; k++) {
        cur_val = tmp[k];
        //std::cerr << "cur val: " << cur_val << std::endl;
        if (cur_val > max_val) {
          max_val = cur_val;
          arg_max = 2*itemp[2*k+1] + k;
        }
      }

      //std::cerr << "minval: " << min_val << std::endl;
      assert(i == nData - (nData % 2));
    }
    
#endif 
    
    for (; i < nData; i++) {
      cur_val = data[i];
      if (cur_val > max_val) {
        max_val = cur_val;
        arg_max = i;
      }
    }
#endif
  }

  inline void find_min_and_argmin(const float_A16* data, const size_t nData, float& min_val, size_t& arg_min)
  {
    min_val = MAX_FLOAT;
    arg_min = MAX_UINT;

    assertAligned16(data);

#if !defined(USE_SSE) || USE_SSE < 4

    if (nData > 0) {
      const float* ptr = std::min_element(data,data+nData);
      min_val = *ptr;
      arg_min = ptr - data;
    }

    // for (i=0; i < nData; i++) {
    //   float cur_val = data[i];

    //   if (cur_val < min_val) {
    //     min_val = cur_val;
    //     arg_min = i;
    //   }
    // }
#else    

    size_t i = 0;
    float cur_val;
    
#if USE_SSE >= 5

    //use AVX  - align16 is no good, need align32 for aligned moves

    if (nData >= 16) {
      const float val = MAX_FLOAT;
      //const uint one = 1;
      const float inc = reinterpret<const uint, const float>(1); //*reinterpret_cast<const float*>(&one);

      asm __volatile__ ("vbroadcastss %[tmp], %%ymm6 \n\t" //ymm6 is min register
                        "vxorps %%ymm5, %%ymm5, %%ymm5 \n\t" //sets ymm5 (= argmin) to zero
                        "vbroadcastss %[itemp], %%ymm4 \n\t" //ymm4 is increment register
                        "vxorps %%ymm3, %%ymm3, %%ymm3 \n\t" //sets ymm3 (= current set index) to zero
                        : : [tmp] "m" (val), [itemp] "m" (inc) : "ymm3", "ymm4", "ymm5", "ymm6");

      for (; (i+8) <= nData; i += 8) {

        asm __volatile__ ("vmovups %[fptr], %%ymm7 \n\t"
                          "vcmpltps %%ymm6, %%ymm7, %%ymm0 \n\t" //destination is last
                          "vblendvps %%ymm0, %%ymm7, %%ymm6, %%ymm6 \n\t" //destination is last
                          "vblendvps %%ymm0, %%ymm3, %%ymm5, %%ymm5 \n\t" //destination is last
                          "vpaddd %%ymm4, %%ymm3, %%ymm3 \n\t" //destination is last
                          : : [fptr] "m" (data[i]) : "ymm0", "ymm3", "ymm5", "ymm6", "ymm7");
      }

      float tmp[8];
      uint itemp[8];

      asm __volatile__ ("vmovups %%ymm6, %[tmp] \n\t"
                        "vmovups %%ymm5, %[itemp]"
                        : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "memory");

      for (uint k=0; k < 8; k++) {
        cur_val = tmp[k];
        //std::cerr << "cur val: " << cur_val << std::endl;
        if (cur_val < min_val) {
          min_val = cur_val;
          arg_min = (itemp[k] << 3) + k; //8*itemp[k] + k;
        }
      }
      assert(i == nData - (nData % 8));
    }

#else
    //blendvps is part of SSE4

    if (nData >= 8) {

      assert(nData <= 17179869183);

      float tmp[4] = {MAX_FLOAT,MAX_FLOAT,MAX_FLOAT,MAX_FLOAT}; //reused as output, static not useful
      uint itemp[4] = {1,1,1,1}; //reused as output, static not useful
      assert(sizeof(uint) == 4);

      asm __volatile__ ("movups %[tmp], %%xmm6 \n\t"
                        "xorps %%xmm5, %%xmm5 \n\t" //sets xmm5 (= argmin) to zero
                        "movups %[itemp], %%xmm4 \n\t"
                        "xorps %%xmm3, %%xmm3 \n\t" //contains candidate argmin
                        : : [tmp] "m" (tmp[0]), [itemp] "m" (itemp[0]) :  "xmm3", "xmm4", "xmm5", "xmm6");

      for (; (i+4) <= nData; i += 4) {

        asm __volatile__ ("movaps %[fptr], %%xmm7 \n\t"
                          "movaps %%xmm7, %%xmm0 \n\t"
                          "cmpltps %%xmm6, %%xmm0 \n\t"
                          "blendvps %%xmm7, %%xmm6 \n\t"
                          "blendvps %%xmm3, %%xmm5 \n\t"
                          "paddd %%xmm4, %%xmm3 \n\t"
                          : : [fptr] "m" (data[i]) : "xmm0", "xmm3", "xmm5", "xmm6", "xmm7");
      }

      asm __volatile__ ("movups %%xmm6, %[tmp] \n\t"
                        "movups %%xmm5, %[itemp] \n\t"
                        : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "memory");

      //std::cerr << "intermediate minval: " << min_val << std::endl;

      for (uint k=0; k < 4; k++) {
        cur_val = tmp[k];
        //std::cerr << "cur val: " << cur_val << std::endl;
        if (cur_val < min_val) {
          min_val = cur_val;
          arg_min = 4*itemp[k] + k;
        }
      }
      assert(i == nData - (nData % 4));
    }

#endif

    //std::cerr << "minval: " << min_val << std::endl;

    for (; i < nData; i++) {
      cur_val = data[i];
      if (cur_val < min_val) {
        min_val = cur_val;
        arg_min = i;
      }
    }
#endif
  }

  inline void find_min_and_argmin(const double_A16* data, const size_t nData, double& min_val, size_t& arg_min)
  {
    min_val = MAX_DOUBLE;
    arg_min = MAX_UINT;

    assertAligned16(data);

#if !defined(USE_SSE) || USE_SSE < 4

    if (nData > 0) {
      const double* ptr = std::min_element(data,data+nData);
      min_val = *ptr;
      arg_min = ptr - data;
    }

    // for (i=0; i < nData; i++) {
    //   double cur_val = data[i];

    //   if (cur_val < min_val) {
    //     min_val = cur_val;
    //     arg_min = i;
    //   }
    // }
#else    

    size_t i = 0;
    double cur_val;
    
#if USE_SSE >= 5

    //use AVX  - align16 is no good, need align32 for aligned moves


    if (nData >= 8) {
      assert(sizeof(size_t) == 8);

      const double val = MAX_DOUBLE;
      //const size_t one = 1;
      const double inc = reinterpret<const size_t, const double>(1); //*reinterpret_cast<const double*>(&one);

      asm __volatile__ ("vbroadcastsd %[tmp], %%ymm6 \n\t" //ymm6 is min register
                        "vxorpd %%ymm5, %%ymm5, %%ymm5 \n\t" //sets ymm5 (= argmin) to zero
                        "vbroadcastsd %[itemp], %%ymm4 \n\t" //ymm4 is increment register
                        "vxorpd %%ymm3, %%ymm3, %%ymm3 \n\t" //sets ymm3 (= current set index) to zero
                        : : [tmp] "m" (val), [itemp] "m" (inc) : "ymm3", "ymm4", "ymm5", "ymm6");

      for (; (i+4) <= nData; i += 4) {

        asm __volatile__ ("vmovupd %[dptr], %%ymm7 \n\t"
                          "vcmpltpd %%ymm6, %%ymm7, %%ymm0 \n\t" //destination is last
                          "vblendvpd %%ymm0, %%ymm7, %%ymm6, %%ymm6 \n\t" //destination is last
                          "vblendvpd %%ymm0, %%ymm3, %%ymm5, %%ymm5 \n\t" //destination is last
                          "vpaddd %%ymm4, %%ymm3, %%ymm3 \n\t" //destination is last
                          : : [dptr] "m" (data[i]) : "ymm0", "ymm3", "ymm5", "ymm6", "ymm7");
      }

      double tmp[4];
      size_t itemp[4];

      asm __volatile__ ("vmovupd %%ymm6, %[tmp] \n\t"
                        "vmovupd %%ymm5, %[itemp]"
                        : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "memory");

      for (uint k=0; k < 4; k++) {
        cur_val = tmp[k];
        //std::cerr << "cur val: " << cur_val << std::endl;
        if (cur_val < min_val) {
          min_val = cur_val;
          arg_min = (itemp[k] << 2) + k; //4*itemp[k] + k;
        }
      }

      assert(i == nData - (nData % 4));
    }

#else

    if (nData >= 4) {

      assert(nData < 8589934592);

      double tmp[2] = {MAX_DOUBLE,MAX_DOUBLE}; //reused as output, static not useful
      uint itemp[4] = {0,1,0,1}; //reused as output, static not useful
      assert(sizeof(uint) == 4);

      //note: broadcast is better implemented by movddup (SSE3). But this code is no longer improved
      asm __volatile__ ("movupd %[tmp], %%xmm6 \n\t"
                        "xorpd %%xmm5, %%xmm5 \n\t" //sets xmm5 (= argmin) to zero
                        "movupd %[itemp], %%xmm4 \n\t"
                        "xorpd %%xmm3, %%xmm3 \n\t" //contains candidate argmin
                        : : [tmp] "m" (tmp[0]), [itemp] "m" (itemp[0]) :  "xmm3", "xmm4", "xmm5", "xmm6");

      for (; (i+2) <= nData; i += 2) {

        asm __volatile__ ("movapd %[dptr], %%xmm7 \n\t"
                          "movapd %%xmm7, %%xmm0 \n\t"
                          "cmpltpd %%xmm6, %%xmm0 \n\t"
                          "blendvpd %%xmm7, %%xmm6 \n\t"
                          "blendvpd %%xmm3, %%xmm5 \n\t"
                          "paddd %%xmm4, %%xmm3 \n\t"
                          : : [dptr] "m" (data[i]) : "xmm0", "xmm3", "xmm5", "xmm6", "xmm7");
      }

      asm __volatile__ ("movupd %%xmm6, %[tmp] \n\t"
                        "movupd %%xmm5, %[itemp] \n\t"
                        : [tmp] "=m" (tmp[0]), [itemp] "=m" (itemp[0]) : : "memory");

      assert(itemp[0] == 0);
      assert(itemp[2] == 0);

      for (uint k=0; k < 2; k++) {
        cur_val = tmp[k];
        //std::cerr << "cur val: " << cur_val << std::endl;
        if (cur_val < min_val) {
          min_val = cur_val;
          arg_min = 2*itemp[2*k+1] + k;
        }
      }
      assert(i == nData - (nData%2));
    }
#endif

    //std::cerr << "minval: " << min_val << std::endl;

    for (; i < nData; i++) {
      cur_val = data[i];
      if (cur_val < min_val) {
        min_val = cur_val;
        arg_min = i;
      }
    }
#endif
  }

  /******************** array mul *********************/

  inline void mul_array(float_A16* data, const size_t nData, const float constant)
  {
    assertAligned16(data);

    size_t i = 0;
#if !defined(USE_SSE) || USE_SSE < 2
    for (; i < nData; i++) { //g++ uses packed avx mul, but after checking alignment
      data[i] *= constant;
    }
#elif USE_SSE >= 5

    // AVX  - align16 is no good, need align32 for aligned moves

    asm __volatile__ ("vbroadcastss %[tmp], %%ymm7 \n\t"
                      : : [tmp] "m" (constant) : "ymm7");

    //vmulps can take an unaligned mem arg!
    for (; i+8 <= nData; i+=8) {
      asm volatile ("vmovups %[fptr], %%ymm6 \n\t"
                    "vmulps %%ymm7, %%ymm6, %%ymm6 \n\t"
                    "vmovups %%ymm6, %[fptr] \n\t"
                    : [fptr] "+m" (data[i]) : : "ymm6", "memory");
    }

    for (; i < nData; i++) 
      data[i] *= constant;
#else
    float temp[4];
    for (; i < 4; i++)
      temp[i] = constant;
    asm volatile ("movaps %[temp], %%xmm7" : : [temp] "m" (temp[0]) : "xmm7" );

    i = 0;
    for (; i+4 <= nData; i+=4) {
      asm volatile ("movaps %[fptr], %%xmm6 \n\t"
                    "mulps %%xmm7, %%xmm6 \n\t"
                    "movaps %%xmm6, %[fptr] \n\t"
                    : [fptr] "+m" (data[i]) : : "xmm6", "memory");
    }

    for (; i < nData; i++) 
      data[i] *= constant;
#endif
  }

  inline void mul_array(double_A16* data, const size_t nData, const double constant)
  {
    size_t i = 0;
#if !defined(USE_SSE) || USE_SSE < 2
    for (; i < nData; i++) {
      data[i] *= constant;
    }
#else
#if USE_SSE >= 5

    // AVX  - align16 is no good, need align32 for aligned moves

    asm __volatile__ ("vbroadcastsd %[tmp], %%ymm7 \n\t"
                      : : [tmp] "m" (constant) : "ymm7");

    //vmulpd can take an unaligned mem arg!
    for (; i+4 <= nData; i+=4) {

      asm volatile ("vmovupd %[dptr], %%ymm6 \n\t"
                    "vmulpd %%ymm7, %%ymm6, %%ymm6 \n\t"
                    "vmovupd %%ymm6, %[dptr] \n\t"
                    : [dptr] "+m" (data[i]) : : "ymm6", "memory");
    }

#else
    double temp[2] = {constant, constant};
    
    //note: broadcast is better implemented by movddup (SSE3). But this code is no longer improved
    asm volatile ("movupd %[temp], %%xmm7" : : [temp] "m" (temp[0]) : "xmm7" );

    for (; i+2 <= nData; i+=2) {

      asm volatile ("movapd %[dptr], %%xmm6 \n\t"
                    "mulpd %%xmm7, %%xmm6 \n\t"
                    "movupd %%xmm6, %[dptr] \n\t"
                    : [dptr] "+m" (data[i]) : : "xmm6", "memory");
    }
#endif

    for (; i < nData; i++)
      data[i] *= constant;
#endif
  }

  /******** array additions with multiplications *************/

  //performs data[i] -= factor*data2[i] for each i
  //this is a frequent operation in the conjugate gradient algorithm
  inline void array_subtract_multiple(double_A16* attr_restrict data, const size_t nData, double factor,
                                      const double_A16* attr_restrict data2)
  {
    assertAligned16(data);

    size_t i = 0;
#if !defined(USE_SSE) || USE_SSE < 2
    for (; i < nData; i++)
      data[i] -= factor * data2[i];
#else    
#if USE_SSE >= 5

    // AVX  - align16 is no good, need align32 for aligned moves

    asm __volatile__ ("vbroadcastsd %[tmp], %%ymm7 \n\t"
                      : : [tmp] "m" (factor) : "ymm7");

    for (; i+4 <= nData; i+=4) {

      //TODO: check for usage of VFNADD

      //vmulpd can take an unaligned mem arg, vfnadd too!
      asm volatile ("vmovupd %[cdptr], %%ymm6 \n\t"
                    "vmulpd %%ymm7, %%ymm6, %%ymm6 \n\t" //destination goes last
                    "vmovupd %[dptr], %%ymm5 \n\t"
                    "vsubpd %%ymm6, %%ymm5, %%ymm5 \n\t" //destination goes last
                    "vmovupd %%ymm5, %[dptr] \n\t"
                    : [dptr] "+m" (data[i]) : [cdptr] "m" (data2[i]) : "ymm5", "ymm6", "memory");
    }
#else
    double temp[2] = {factor, factor};
    
    //note: broadcast is better implemented by movddup (SSE3). But this code is no longer improved
    asm volatile ("movupd %[temp], %%xmm7" : : [temp] "m" (temp[0]) : "xmm7" );
    for (; i+2 <= nData; i+=2) {

      asm volatile ("movapd %[cdptr], %%xmm6 \n\t"
                    "mulpd %%xmm7, %%xmm6 \n\t"
                    "movapd %[dptr], %%xmm5 \n\t"
                    "subpd %%xmm6, %%xmm5 \n\t"
                    "movapd %%xmm5, %[dptr] \n\t"
                    : [dptr] "+m" (data[i]) : [cdptr] "m" (data2[i]) : "xmm5", "xmm6", "memory");
    }
#endif

    for (; i < nData; i++)
      data[i] -= factor * data2[i];
#endif
  }

  inline void array_add_multiple(double_A16* attr_restrict data, const size_t nData, double factor,
                                 const double_A16* attr_restrict data2)
  {
    array_subtract_multiple(data, nData, -factor, data2);
  }

  //NOTE: despite attr_restrict, you can safely pass the same for dest and src1 or src2
  inline void go_in_neg_direction(double_A16* attr_restrict dest, const size_t nData, const double_A16* attr_restrict src1,
                                  const double_A16* attr_restrict src2, double alpha)
  {
    assertAligned16(dest);

    size_t i = 0;
#if !defined(USE_SSE) || USE_SSE < 5
    for (; i < nData; i++)
      dest[i] = src1[i] - alpha * src2[i];
#else

    // AVX  - align16 is no good, need align32 for aligned moves

    asm __volatile__ ("vbroadcastsd %[w1], %%ymm0 \n\t" //ymm0 = w1
                      : : [w1] "m" (alpha) : "ymm0");

    //vmulpd can take an unaligned mem arg, vfnmadd too!
    //vsubpd can take an unaligned mem arg!
    for (; i+4 <= nData; i+= 4) {

      asm volatile ("vmovupd %[s2_ptr], %%ymm3 \n\t"
#if USE_SSE < 6
                    "vmulpd %%ymm0, %%ymm3, %%ymm3 \n\t" //destination goes last
                    "vmovupd %[s1_ptr], %%ymm2 \n\t"
                    "vsubpd %%ymm3, %%ymm2, %%ymm2 \n\t" //destination goes last
#else
                    "vmovupd %[s1_ptr], %%ymm2 \n\t"
                    "vfnmadd231pd %%ymm0, %%ymm3, %%ymm2 \n\t" //destination goes last
#endif
                    "vmovupd %%ymm2, %[dest]"
                    : [dest] "+m" (dest[i]) : [s1_ptr] "m" (src1[i]), [s2_ptr] "m" (src2[i]) : "ymm2", "ymm3", "memory");
    }

    for (; i < nData; i++)
      dest[i] = src1[i] - alpha * src2[i];
#endif
  }

  //NOTE: despite attr_restrict, you can safely pass the same for dest and src1 or src2
  inline void assign_weighted_combination(double_A16* attr_restrict dest, const size_t nData, double w1, const double_A16* attr_restrict src1,
                                          double w2, const double_A16* attr_restrict src2)
  {
    size_t i = 0;
#if !defined(USE_SSE) || USE_SSE < 5
    for (; i < nData; i++)
      dest[i] = w1 * src1[i] + w2 * src2[i];
#else
    //use AVX

    asm __volatile__ ("vbroadcastsd %[w1], %%ymm0 \n\t" //ymm0 = w1
                      "vbroadcastsd %[w2], %%ymm1 \n\t" //ymm1 = w2
                      : : [w1] "m" (w1), [w2] "m" (w2) : "ymm0", "ymm1");

    //vmulpd can take an unaligned mem arg, vfnmadd too!
    //vaddpd can take an unaligned mem arg!
    for (; i+4 <= nData; i+= 4) {

      asm volatile ("vmovupd %[s1_ptr], %%ymm2 \n\t"
                    "vmulpd %%ymm0, %%ymm2, %%ymm2 \n\t" //destination goes last
                    "vmovupd %[s2_ptr], %%ymm3 \n\t"
#if USE_SSE < 6
                    "vmulpd %%ymm1, %%ymm3, %%ymm3 \n\t" //destination goes last
                    "vaddpd %%ymm3, %%ymm2, %%ymm2 \n\t" //destination goes last
#else
                    "vfmadd231pd %%ymm3, %%ymm1, %%ymm2 \n\t" //destination goes last
#endif
                    "vmovupd %%ymm2, %[dest]"
                    : [dest] "+m" (dest[i]) : [s1_ptr] "m" (src1[i]), [s2_ptr] "m" (src2[i]) : "ymm2", "ymm3", "memory");
    }

    for (; i < nData; i++)
      dest[i] = w1 * src1[i] + w2 * src2[i];
#endif
  }

  /***************** dot product  *****************/  
  
  template<typename T>
  inline T dotprod(const T* data1, const T* data2, const size_t size) 
  {
    return std::inner_product(data1, data1+size, data2, (T) 0);
  }

  template<>
  inline double dotprod(const double* data1, const double* data2, const size_t size) 
  {
#if !defined(USE_SSE) || USE_SSE < 5          
    return std::inner_product((double_A16*) data1, (double_A16*) data1+size, data2, 0.0);
#else     
    //checked: g++ does not use dppd. It uses 256 bit instead, but the running times are the same

    //NOTE: unlike vpps, vppd is not available for 256 bit, not even in AVX-512

    asm __volatile__ ("vxorpd %%xmm12, %%xmm12, %%xmm12 \n\t" : : : "xmm12"); 

    double result = 0.0;

    size_t i = 0;
    for (; i + 2 <= size; i += 2) 
    {
      //std::cerr << "i: " << i << std::endl;  
        
      //vdppd can take a mem arg!
      asm __volatile__ ("vmovupd %[d1], %%xmm10 \n\t"
                        "vmovupd %[d2], %%xmm11 \n\t"
                        "vdppd $49, %%xmm10, %%xmm11, %%xmm13 \n\t" //include all, write in first (hence second is set to 0)
                        "vaddsd %%xmm13, %%xmm12, %%xmm12 \n\t"
                        : : [d1] "m" (data1[i]), [d2] "m" (data2[i]) : "xmm10", "xmm11", "xmm12", "xmm13");
                        
      //std::cerr << "state after add: " << temp[0] << "," << temp[1] << std::endl;
    }

    asm __volatile__ ("vmovlpd %%xmm12, %0 \n\t" : "+m" (result) : : );

    for (; i < size; i++)
      result += data1[i] * data2[i];
    
    return result;
#endif
  }

  /******************** binary search *********************/

  //binary search, returns MAX_UINT if key is not found, otherwise the position in the vector
  template<typename T>
  inline size_t binsearch(const T* data, const T key, const size_t nData)
  {
    if (nData == 0 || key < data[0] || key > data[nData-1])
      return MAX_UINT;

    size_t lower = 0;
    size_t upper = nData-1;
    if (data[lower] == key)
      return lower;
    if (data[upper] == key)
      return upper;

    while (lower+1 < upper) {
      assert(data[lower] < key);
      assert(data[upper] > key);

      const size_t middle = (lower+upper) >> 1;  // (lower+upper)/2;
      assert(middle > lower && middle < upper);
      if (data[middle] == key)
        return middle;
      else if (data[middle] < key)
        lower = middle;
      else
        upper = middle;
    }

    return MAX_UINT;
  }

  template<typename T>
  inline size_t binsearch_insertpos(const T* data, const T key, const size_t nData)
  {
    if (nData == 0 || key <= data[0])
      return 0;

    if (key > data[nData-1])
      return nData;

    size_t lower = 0;
    size_t upper = nData-1;
    if (data[upper] == key)
      return upper;

    while (lower+1 < upper) {
      assert(data[lower] < key);
      assert(data[upper] > key);

      const size_t middle = (lower+upper) >> 1;  // (lower+upper)/2;
      assert(middle > lower && middle < upper);
      if (data[middle] == key)
        return middle;
      else if (data[middle] < key)
        lower = middle;
      else
        upper = middle;
    }

    assert(lower+1 == upper);
    return upper;
  }
}

#endif