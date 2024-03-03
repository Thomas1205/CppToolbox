/**************** written by Thomas Schoenemann, January 2024 *******************/
/**************** vector class without virtual methods, saving the pointer to vtable. Useful when you have many vectors **************/


#ifndef PLAIN_VECTOR_HH
#define PLAIN_VECTOR_HH

#include "plain_storage1D.hh"
#include "makros.hh"

template<typename T, typename ST = uchar> //if you want to save the vtable entry, most likely your vector will fit into 255 bytes
class PlainVector : public PlainStorage1D<T,ST> {
public:	

  using Base = PlainStorage1D<T,ST>;	

  //according to https://gcc.gnu.org/onlinedocs/gcc-7.2.0/gcc/Common-Type-Attributes.html#Common-Type-Attributes , alignment has to be expressed like this:
  typedef T T_A16 ALIGNED16;

  explicit PlainVector() noexcept : PlainStorage1D<T,ST>() {}
 
  explicit PlainVector(ST size) : PlainStorage1D<T,ST>(size) {}
 
  explicit PlainVector(ST size, T default_value) : PlainStorage1D<T,ST>(size,default_value) {}
 
  PlainVector(const std::initializer_list<T>& list) : PlainStorage1D<T,ST>(list) {}

  PlainVector(const PlainVector<T,ST>& toCopy) : PlainStorage1D<T,ST>(toCopy) {}
 
  PlainVector(PlainVector<T,ST>&& toTake) noexcept : PlainStorage1D<T,ST>(toTake) {}
  
  PlainVector<T,ST>& operator=(const PlainVector<T,ST>& toCopy) 
  {
	delete[] Base::data_;
	if (toCopy.data_ == 0)
	{
	  Base::data_ = 0;
	  Base::size_ = 0;
	  assert(toCopy.size_ == 0);
	}
	else {
  	  Base::size_ = toCopy.size_;
	  Base::data_ = new T[Base::size_];
      Makros::unified_assign(Base::data_, toCopy.direct_access(), Base::size_);	  
	}
	return *this;  
  }

  PlainVector<T,ST>& operator=(PlainVector<T,ST>&& toTake) noexcept
  {
    delete[] Base::data_;
    Base::data_ = toTake.data_;
    Base::size_ = toTake.size_;
    toTake.data_ = 0;

    return *this;
  }  

  inline T sum() const noexcept
  {
     const ST size = Base::size_;
    const T_A16* data = Base::data_;

    assertAligned16(data);

    return std::accumulate(data,data+size,(T)0); 
  }
  
  /*** maximal element ***/
  T max() const noexcept {
    const ST size = Base::size_;

    if (size > 0) {
      // g++ 4.8.5 uses avx instructions for this
      const T_A16* data = Base::data_;
      assertAligned16(data);

      return *std::max_element(data, data + size);
    }
    else
      return std::numeric_limits<T>::min();
  }

  /*** minimal element ***/
  T min() const noexcept {
    const ST size = Base::size_;

    if (size > 0) {
      // g++ 4.8.5 uses avx instructions for this
      const T_A16* data = Base::data_;
      assertAligned16(data);

      return *std::min_element(data, data + size);
    }
    else
      return std::numeric_limits<T>::max();
  }
  
  void operator+=(const PlainVector<T,ST>& v) noexcept
  {
    const ST size = Base::size_;

#ifndef DONT_CHECK_VECTOR_ARITHMETIC
    if (v.size() != size) {
      INTERNAL_ERROR << "   cannot add PlainVector to PlainVector: "
                     << "   sizes " << v.size() << " and " << this->size() << " mismatch. Exiting..."
                     << std::endl;
      exit(0);
    }
#endif

    ST i;

    // for (i=0; i < size; i++)
    //   Base::data_[i] += v.direct_access(i);

    T_A16* attr_restrict dptr = Base::data_;
    const T_A16* attr_restrict vptr = v.direct_access();

    assertAligned16(dptr);
    assertAligned16(vptr);

    for (i=0; i < size; i++)
      dptr[i] += vptr[i];	  
  }

  void operator-=(const PlainVector<T,ST>& v) noexcept
  {
    const ST size = Base::size_;

#ifndef DONT_CHECK_VECTOR_ARITHMETIC
    if (v.size() != size) {
      INTERNAL_ERROR << "   cannot add PlainVector to PlainVector: "
                     << "   sizes " << v.size() << " and " << this->size() << " mismatch. Exiting..."
                     << std::endl;
      exit(0);
    }
#endif

    ST i;
    // for (i=0; i < size; i++)
    //   Base::data_[i] -= v.direct_access(i);

    T_A16* attr_restrict dptr = Base::data_;
    const T_A16* attr_restrict vptr = v.direct_access();

    for (i=0; i < size; i++)
      dptr[i] -= vptr[i];	  
  }

  void operator*=(const T constant) noexcept
  {
    const ST size = Base::size_;
    T_A16* data = Base::data_;

    assertAligned16(data);

    ST i;
    for (i=0; i < size; i++)
      data[i] *= constant;	  
  }
  
};

template<typename T, typename ST>
T operator%(const PlainVector<T,ST>& v1, const PlainVector<T,ST>& v2) {
  const ST size = v1.size();	

#ifndef DONT_CHECK_VECTOR_ARITHMETIC
    if (size != v2.size()) {
      INTERNAL_ERROR << "     cannot compute scalar product of PlainVectors: \"" << std::endl
                     << "      sizes " << v1.size() << " and " << v2.size() << " mismatch. exiting."
                     << std::endl;
      exit(1);
    }
#endif
	
  return Routines::dotprod(v1.direct_access(),v2.direct_access(),size);	
}

#endif