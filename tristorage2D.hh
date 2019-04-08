/*-*-c++-*-*/ 
/*** written by Thomas Schoenemann as a private person without employment, March 2015 ***/
/*** if you desire the checked version, make sure your compiler defines the option SAFE_MODE on the command line ***/

#ifndef TRISTORAGE2D_HH
#define TRISTORAGE2D_HH

#include "makros.hh"

//two-dimensional container class for objects of any type T, where only have the pattern can be filled/defined
//you can use this either for a triangular pattern or for a symmetric one 
//(as long as we use only the storage it does not matter much, except when printing. for matrix usage it does matter, e.g. in the matrix-vector product)
//(neither mathematical nor streaming operations need to be defined on T)
template<typename T, typename ST = size_t>
class TriStorage2D {
public:

  //default constructor
  TriStorage2D();

  TriStorage2D(ST dim);

  TriStorage2D(ST dim, T default_value);

  //copy constructor
  TriStorage2D(const TriStorage2D<T,ST>& toCopy);

  ~TriStorage2D();

  virtual const std::string& name() const;

  //existing elements are preserved, new ones are uninitialized
  void resize(ST newDim);

  //existing elements are preserved, new ones are filled with the second argument
  void resize(ST newDim, const T fill_value);

  //all elements are uninitialized after this operation
  void resize_dirty(ST newDim);

  void set_constant(const T new_constant);

  //access on an element (handling is symmetric, i.e. accessing (x,y) is equivalent to accessing (y,x) )
  OPTINLINE const T& operator()(ST x, ST y) const;

  OPTINLINE T& operator()(ST x, ST y);

  void operator=(const TriStorage2D<T,ST>& toCopy);

  inline T* direct_access();

  inline const T* direct_access() const;

  inline T& direct_access(ST i);

  inline T direct_access(ST i) const;

  inline ST dim() const;
  
  inline ST size() const;

protected:

  T* data_;
  ST dim_;
  ST size_;
  static const std::string tristor2D_name_;
};


template<typename T, typename ST=size_t>
class NamedTriStorage2D : public TriStorage2D<T,ST> {
public:

  NamedTriStorage2D();

  NamedTriStorage2D(std::string name);
  
  NamedTriStorage2D(ST dim, std::string name);
  
  NamedTriStorage2D(ST dim, T default_value, std::string name);

  virtual const std::string& name() const;

  inline void operator=(const TriStorage2D<T,ST>& toCopy);

  //NOTE: the name is NOT copied
  inline void operator=(const NamedTriStorage2D<T,ST>& toCopy);

protected:
  std::string name_;
};


template<typename T, typename ST>
bool operator==(const TriStorage2D<T,ST>& v1, const TriStorage2D<T,ST>& v2);

template<typename T, typename ST>
bool operator!=(const TriStorage2D<T,ST>& v1, const TriStorage2D<T,ST>& v2);

namespace Makros {

  template<typename T, typename ST>
  class Typename<TriStorage2D<T,ST> > {
  public:

    std::string name() const {

      return "TriStorage2D<" + Typename<T>() + "," + Typename<ST>() + "> ";
    }
  };

  template<typename T>
  class Typename<TriStorage2D<T> > {
  public:

    std::string name() const {

      return "TriStorage2D<" + Typename<T>() + "> ";
    }
  };

  template<typename T, typename ST>
  class Typename<NamedTriStorage2D<T,ST> > {
  public:

    std::string name() const {

      return "NamedTriStorage2D<" + Typename<T>() + "," + Typename<ST>() + "> ";
    }
  };

  template<typename T>
  class Typename<NamedTriStorage2D<T> > {
  public:

    std::string name() const {

      return "NamedTriStorage2D<" + Typename<T>() + "> ";
    }
  };

}


/**************************** implementation **************************************/

template<typename T, typename ST>
/*static*/ const std::string TriStorage2D<T,ST>::tristor2D_name_ = "unnamed TriStorage2D";

//constructors
template<typename T, typename ST>
TriStorage2D<T,ST>::TriStorage2D() : data_(0), dim_(0), size_(0) {}

template<typename T, typename ST>
TriStorage2D<T,ST>::TriStorage2D(ST dim) : dim_(dim) {

  size_ = dim_*(dim_+1) / 2;
  data_ = new T[size_];
}

template<typename T, typename ST>
TriStorage2D<T,ST>::TriStorage2D(ST dim, T default_value) : dim_(dim) {

  size_ = dim_*(dim_+1) / 2;
  data_ = new T[size_];

  std::fill(data_, data_+size_, default_value); //fill and fill_n are of equal speed
}

//copy constructor
template<typename T, typename ST>
TriStorage2D<T,ST>::TriStorage2D(const TriStorage2D<T,ST>& toCopy) {

  dim_ = toCopy.dim();
  size_ = toCopy.size();

  if (size_ == 0)
    data_ = 0;
  else {
    data_ = new T[size_];

    for (ST i = 0; i < size_; i++)
      data_[i] = toCopy.direct_access(i);
  }
}

//destructor
template <typename T, typename ST>
TriStorage2D<T,ST>::~TriStorage2D() {
  if (data_ != 0)
    delete[] data_;
}

template<typename T, typename ST>
void TriStorage2D<T,ST>::operator=(const TriStorage2D<T,ST>& toCopy) {

  dim_ = toCopy.dim();
  size_ = toCopy.size();

  //room for improvement here: could check if the dimensions already match, then we can reuse data_

  if (data_ != 0)
    delete[] data_;

  if (size_ == 0)
    data_ = 0;
  else {
    data_ = new T[size_];

    for (ST i = 0; i < size_; i++)
      data_[i] = toCopy.direct_access(i);
  }
}

template <typename T, typename ST>
/*virtual*/ const std::string& TriStorage2D<T,ST>::name() const {
  return tristor2D_name_;
}

template <typename T, typename ST>
void TriStorage2D<T,ST>::set_constant(const T new_constant) {
  
  std::fill_n(data_,size_,new_constant);

  // for (ST i=0; i < size_; i++)
  //   data_[i] = new_constant;
}


template<typename T, typename ST>
inline T* TriStorage2D<T,ST>::direct_access() { return data_; }

template<typename T, typename ST>
inline const T* TriStorage2D<T,ST>::direct_access() const { return data_; }

template<typename T, typename ST>
inline T& TriStorage2D<T,ST>::direct_access(ST i) { return data_[i]; }

template<typename T, typename ST>
inline T TriStorage2D<T,ST>::direct_access(ST i) const {
  return data_[i];
}

template<typename T, typename ST>
inline ST TriStorage2D<T,ST>::dim() const { return dim_; }

template<typename T, typename ST>
inline ST TriStorage2D<T,ST>::size() const { return size_; }


template <typename T, typename ST>
OPTINLINE const T& TriStorage2D<T,ST>::operator()(ST x, ST y) const {
#ifdef SAFE_MODE
  if (x >= dim_ || y >= dim_) {
    INTERNAL_ERROR << "    access on element(" << x << "," << y 
                   << ") exceeds storage dimensions of (" << dim_ << "," << dim_ << ")" << std::endl;
    std::cerr << "      in TriStorage2D \"" << this->name() << "\" of type " 
              << Makros::Typename<T>()
              << ". Exiting." << std::endl;
    exit(1);
  }
#endif
  if (x > y)
    std::swap(x,y);

  return data_[(y*(y+1))/2+x];
}

template <typename T, typename ST>
OPTINLINE T& TriStorage2D<T,ST>::operator()(ST x, ST y)  {
#ifdef SAFE_MODE
  if (x >= dim_ || y >= dim_) {
    INTERNAL_ERROR << "    access on element(" << x << "," << y 
                   << ") exceeds storage dimensions of (" << dim_ << "," << dim_ << ")" << std::endl;
    std::cerr << "      in TriStorage2D \"" << this->name() << "\" of type " 
              << Makros::Typename<T>()
              << ". Exiting." << std::endl;
    exit(1);
  }
#endif

  //std::cerr << " access(" << x << "," << y << ")"; // << std::endl;
  
  if (x > y) 
    std::swap(x,y);

  return data_[(y*(y+1))/2+x];
}


//existing elements are preserved, new ones are uninitialized
template<typename T, typename ST>
void TriStorage2D<T,ST>::resize(ST newDim) {

  if (newDim != dim_) {

    ST new_size = newDim*(newDim+1) / 2;
    T* new_data = new T[new_size];

    //copy existing values
    for (ST i=0; i < std::min(size_,new_size); i++)
      new_data[i] = data_[i];

    if (data_ != 0)
      delete[] data_;

    data_ = new_data;
    size_ = new_size;
    dim_ = newDim;
  }
}


//existing elements are preserved, new ones are filled with the second argument
template<typename T, typename ST>
void TriStorage2D<T,ST>::resize(ST newDim, const T fill_value) {

  if (newDim != dim_) {

    ST new_size = newDim*(newDim+1) / 2;
    T* new_data = new T[new_size];

    //copy existing values
    for (ST i=0; i < std::min(size_,new_size); i++)
      new_data[i] = data_[i];

    //fill new values
    if (new_size > size_)
      std::fill_n(new_data+size_,new_size-size_,fill_value);
    // for (ST i=std::min(size_,new_size); i < new_size; i++)
    //   new_data[i] = fill_value;

    if (data_ != 0)
      delete[] data_;

    data_ = new_data;
    size_ = new_size;
    dim_ = newDim;
  }
}

template<typename T, typename ST>
void TriStorage2D<T,ST>::resize_dirty(ST newDim) {

  if (newDim != dim_) {
    if (data_ != 0) {
      delete[] data_;
    }
  
    dim_ = newDim;
    size_ = dim_*(dim_+1) / 2;
    data_ = new T[size_];
  }
}


template<typename T, typename ST>
bool operator==(const TriStorage2D<T,ST>& v1, const TriStorage2D<T,ST>& v2) {

  if (v1.size() != v2.size())
    return false;
  
  uint i;
  for (i=0; i < v1.size(); i++) {
    if (v1.direct_access(i) != v2.direct_access(i)) 
      return false;
  }

  return true;
}

template<typename T, typename ST>
bool operator!=(const TriStorage2D<T,ST>& v1, const TriStorage2D<T,ST>& v2) {

  if (v1.size() != v2.size())
    return true;
  
  uint i;
  for (i=0; i < v1.size(); i++) {
    if (v1.direct_access(i) != v2.direct_access(i)) 
      return true;
  }

  return false;
}

/***** implementation of NamedTriStorage2D ********/

template<typename T, typename ST>
NamedTriStorage2D<T,ST>::NamedTriStorage2D() {}

template<typename T, typename ST>
NamedTriStorage2D<T,ST>::NamedTriStorage2D(std::string name) : name_(name) {}
  
template<typename T, typename ST>
NamedTriStorage2D<T,ST>::NamedTriStorage2D(ST dim, std::string name) : TriStorage2D<T,ST>(dim), name_(name) {}
  
template<typename T, typename ST>
NamedTriStorage2D<T,ST>::NamedTriStorage2D(ST dim, T default_value, std::string name) : TriStorage2D<T,ST>(dim,default_value), name_(name) {}

template<typename T, typename ST>
/*virtual*/ const std::string& NamedTriStorage2D<T,ST>::name() const {
  return name_;
}

template<typename T, typename ST>
inline void NamedTriStorage2D<T,ST>::operator=(const TriStorage2D<T,ST>& toCopy) {
  TriStorage2D<T,ST>::operator=(toCopy);
}

//NOTE: the name is NOT copied
template<typename T, typename ST>
inline void NamedTriStorage2D<T,ST>::operator=(const NamedTriStorage2D<T,ST>& toCopy) {
  TriStorage2D<T,ST>::operator=(toCopy);
}


#endif
