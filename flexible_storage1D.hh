/********* Written by Thomas Schoenemann, created from moved code, February 2020 **********/

#ifndef FLEXIBLE_STORAGE1D_HH
#define FLEXIBLE_STORAGE1D_HH

#include "makros.hh"

//this class is meant to replace std::vector with its push_back() functionality.
// It has slightly less functionality, though. E.g. erase() is not available.
template<typename T, typename ST=size_t>
class FlexibleStorage1D {
public:

  FlexibleStorage1D();

  FlexibleStorage1D(ST reserved_size);

  //copy constructor
  FlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy);

  ~FlexibleStorage1D();

  virtual const std::string& name() const;

  inline T& operator[](ST i) const;

  inline ST size() const;

  inline ST reserved_size() const;

  //for compatibility with std::vector, e.g. for use in templates
  inline ST capacity() const;

  void resize(ST size, bool exact_fit = false);

  void fit_exactly();

  void shrink(ST size);

  inline void shrink_by(ST reduction);
  
  inline void grow_by_dirty(ST increase);

  void reserve(ST size);

  //will not free memory
  void clear();

  void set_constant(T val);

  void operator=(const FlexibleStorage1D<T,ST>& toCopy);

  inline T back() const;

  inline T& back();

  ST append(T val);

  //shortcut when you are sure the allocated memory suffices
  inline void append_trusting(T val);

  //for compatibility with std::vector, e.g. for use in templates
  inline void push_back(T val);

  void append(Storage1D<T,ST>& toAppend);

  void append(FlexibleStorage1D<T,ST>& toAppend);

  void erase(ST pos);
  
  void erase_several(ST pos, ST nToErase);

  void insert(ST pos, T val);

  T* direct_access();

  const T* direct_access() const;

  void swap(FlexibleStorage1D<T,ST>& toSwap);

protected:

  T* data_;
  ST size_;
  ST reserved_size_;
  static const std::string flex_stor1D_name_;
};

template<typename T, typename ST>
std::ostream& operator<<(std::ostream& s, const FlexibleStorage1D<T,ST>& v);

template<typename T, typename ST>
bool operator==(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2);

template<typename T, typename ST>
bool operator!=(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2);

template<typename T, typename ST=size_t>
class NamedFlexibleStorage1D : public FlexibleStorage1D<T,ST> {
public:

  NamedFlexibleStorage1D();

  NamedFlexibleStorage1D(const std::string& name);

  NamedFlexibleStorage1D(ST reserved_size, const std::string& name);

  //copy constructors
  NamedFlexibleStorage1D(const NamedFlexibleStorage1D<T,ST>& toCopy);

  NamedFlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy);

  virtual const std::string& name() const;

  //operators
  void operator=(const NamedFlexibleStorage1D<T,ST>& toCopy);

  void operator=(const FlexibleStorage1D<T,ST>& toCopy);

protected:
  std::string name_;
};

/******* implementation of FlexibleStorage1D *********/

template<typename T, typename ST>
/*static*/ const std::string FlexibleStorage1D<T,ST>::flex_stor1D_name_ = "unnamed flexible 1Dstorage";

template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::FlexibleStorage1D() : size_(0)
{
  reserved_size_ = 4;
  data_ = new T[reserved_size_];
}

template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::FlexibleStorage1D(ST reserved_size)  : size_(0), reserved_size_(reserved_size)
{
  data_ = new T[reserved_size_];
}

//copy constructor
template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::FlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy)
{
  size_ = toCopy.size();
  reserved_size_ = toCopy.reserved_size();

  data_ = new T[reserved_size_];

  Makros::unified_assign(data_, toCopy.direct_access(), size_);

  //for (uint k=0; k < toCopy.size(); k++)
  //  data_[k] = toCopy[k];
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::swap(FlexibleStorage1D<T,ST>& toSwap)
{
  std::swap(data_,toSwap.data_);
  std::swap(size_,toSwap.size_);
  std::swap(reserved_size_,toSwap.reserved_size_);
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::operator=(const FlexibleStorage1D<T,ST>& toCopy)
{
  uint new_res = toCopy.reserved_size();
  if (new_res != reserved_size_) {
    reserved_size_ = new_res;

    if (data_ != 0)
      delete[] data_;
    data_ = new T[reserved_size_];
  }

  size_ = toCopy.size();

  Makros::unified_assign(data_, toCopy.direct_access(), size_);

  //for (uint k=0; k < size_; k++)
  //  data_[k] = toCopy[k];
}

template<typename T, typename ST>
/*virtual*/ const std::string& FlexibleStorage1D<T,ST>::name() const
{
  return flex_stor1D_name_;
}

template<typename T, typename ST> 
FlexibleStorage1D<T,ST>::~FlexibleStorage1D()
{
  delete[] data_;
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::set_constant(T val)
{
  for (ST k=0; k < size_; k++)
    data_[k] = val;
}

template<typename T, typename ST>
inline ST FlexibleStorage1D<T,ST>::size() const
{
  return size_;
}

template<typename T, typename ST>
inline ST FlexibleStorage1D<T,ST>::reserved_size() const
{
  return reserved_size_;
}

template<typename T, typename ST>
inline ST FlexibleStorage1D<T,ST>::capacity() const
{
  return reserved_size_;
}

template<typename T, typename ST>
inline T FlexibleStorage1D<T,ST>::back() const
{
  assert(size_ > 0);
  assert(data_ != 0);
  return data_[size_-1];
}

template<typename T, typename ST>
inline T& FlexibleStorage1D<T,ST>::back()
{
  assert(size_ > 0);
  assert(data_ != 0);
  return data_[size_-1];
}

template<typename T, typename ST>
ST FlexibleStorage1D<T,ST>::append(T val)
{
  if (size_ == reserved_size_) {

    reserved_size_ = size_t(1.2 * reserved_size_) + 4;

    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, data_, size_);

    //for (uint k=0; k < size_; k++)
    //  new_data[k] = data_[k];

    delete[] data_;
    data_ = new_data;
  }

  const ST k = size_;
  data_[k] = val;

  size_++;

  return k;
}

//shortcut when you are sure the allocated memory suffices
template<typename T, typename ST>
inline void FlexibleStorage1D<T,ST>::append_trusting(T val)
{
  assert(size_ < reserved_size_);
  data_[size_] = val;
  size_++;
}

template<typename T, typename ST>
inline void FlexibleStorage1D<T,ST>::push_back(T val)
{
  append(val);
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::append(Storage1D<T,ST>& toAppend)
{
  if (reserved_size_ < size_ + toAppend.size()) {

    reserved_size_ = size_ + toAppend.size() + 2;

    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, data_, size_);

    //for (uint k=0; k < size_; k++)
    //  new_data[k] = data_[k];

    delete[] data_;
    data_ = new_data;
  }

  for (ST k=0; k < toAppend.size(); k++) {
    data_[size_] = toAppend[k];
    size_++;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::append(FlexibleStorage1D<T,ST>& toAppend)
{
  if (reserved_size_ < size_ + toAppend.size()) {

    reserved_size_ = size_ + toAppend.size() + 2;

    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, data_, size_);

    //for (uint k=0; k < size_; k++)
    //  new_data[k] = data_[k];

    delete[] data_;
    data_ = new_data;
  }

  for (ST k=0; k < toAppend.size(); k++) {
    data_[size_] = toAppend[k];
    size_++;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::resize(ST size, bool exact_fit)
{
  if (size > reserved_size_ || size < (reserved_size_ / 3) ) {

    reserved_size_ = size;
    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, data_, std::min(size_,size));

    //for (uint k=0; k < std::min(size_,size); k++)
    //  new_data[k] = data_[k];

    delete[] data_;
    data_ = new_data;
  }

  size_ = size;

  if (exact_fit && size_ != reserved_size_) {

    reserved_size_ = size_;
    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, data_, size_);

    //for (uint k=0; k < size_; k++)
    //  new_data[k] = data_[k];

    delete[] data_;
    data_ = new_data;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::fit_exactly()
{
  if (reserved_size_ != size_) {
    
    T* new_data = new T[size_];
    Makros::unified_assign(new_data, data_, size_);
    reserved_size_ = size_;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::shrink(ST size)
{
  assert(size <= size_);
  size_ = size;
}

template<typename T, typename ST>
inline void FlexibleStorage1D<T,ST>::shrink_by(ST reduction)
{
  assert(reduction <= size_);
  size_ -= reduction;
}

template<typename T, typename ST>
inline void FlexibleStorage1D<T,ST>::grow_by_dirty(ST increase)
{
  const ST size = size_ + increase;
 
  if (size > reserved_size_ || size < (reserved_size_ / 3) ) {

    reserved_size_ = size;
    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, data_, std::min(size_,size));

    //for (uint k=0; k < std::min(size_,size); k++)
    //  new_data[k] = data_[k];

    delete[] data_;
    data_ = new_data;
  }

  size_ = size; 
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::reserve(ST size)
{
  if (size > size_ && size != reserved_size_) {

    reserved_size_ = size;
    T* new_data = new T[reserved_size_];

    Makros::unified_assign(new_data, data_, size_);

    delete[] data_;
    data_ = new_data;
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::clear()
{
  size_ = 0;
}

template<typename T, typename ST>
inline T& FlexibleStorage1D<T,ST>::operator[](ST i) const
{

#ifdef SAFE_MODE
  if (i >= size_) {
    INTERNAL_ERROR << "    invalid access on element " << i
                   << " for FlexibleStorage1D " <<  "\"" << this->name() << "\" of type "
                   //<< Makros::Typename<T>()
                   << typeid(T).name()
                   << " with " << size_ << " (valid) elements. exiting." << std::endl;
    print_trace();
    exit(1);
  }
#endif
  return data_[i];
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::erase(ST pos)
{
  if (pos < size_)
  {
    Makros::downshift_array(data_, pos, size_, 1);
    shrink_by(1);
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::erase_several(ST pos, ST nToErase)
{
  if (pos < size_) 
  {
    nToErase = std::min(nToErase, size_ - pos);
    Makros::downshift_array(data_, pos, size_, nToErase);
    shrink_by(nToErase);    
  }
}

template<typename T, typename ST>
void FlexibleStorage1D<T,ST>::insert(ST pos, T val)
{
  if (pos <= size_)
  {
    grow_by_dirty(1);
    Makros::upshift_array(data_, pos, size_, 1);
    data_[pos] = val;
  }
}

template<typename T, typename ST>
T* FlexibleStorage1D<T,ST>::direct_access()
{
  return data_;
}

template<typename T, typename ST>
const T* FlexibleStorage1D<T,ST>::direct_access() const
{
  return data_;
}

template<typename T, typename ST>
std::ostream& operator<<(std::ostream& s, const FlexibleStorage1D<T,ST>& v)
{
  s << "[ ";
  for (int i=0; i < ((int) v.size()) - 1; i++)
    s << v[i] << ",";
  if (v.size() > 0)
    s << v[v.size()-1];
  s << " ]";

  return s;
}

template<typename T, typename ST>
bool operator==(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2)
{
  if (v1.size() != v2.size())
    return false;

  for (ST k=0; k < v1.size(); k++) {
    if (v1[k] != v2[k])
      return false;
  }
  return true;
}

template<typename T, typename ST>
bool operator!=(const FlexibleStorage1D<T,ST>& v1, const FlexibleStorage1D<T,ST>& v2)
{
  if (v1.size() != v2.size())
    return true;

  for (ST k=0; k < v1.size(); k++) {
    if (v1[k] != v2[k])
      return true;
  }
  return false;
}

/***********************************/

template<typename T, typename ST> NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D() : name_("unfs1d") {}

template<typename T, typename ST> NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(const std::string& name) : name_(name)
{
}

template<typename T, typename ST> NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(ST reserved_size, const std::string& name) :
  FlexibleStorage1D<T,ST>(reserved_size), name_(name) {}

//Note: the name is NOT copied
template<typename T, typename ST> NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(const NamedFlexibleStorage1D<T,ST>& toCopy) :
  FlexibleStorage1D<T,ST>(toCopy), name_("unfs1d")
{
}

template<typename T, typename ST> NamedFlexibleStorage1D<T,ST>::NamedFlexibleStorage1D(const FlexibleStorage1D<T,ST>& toCopy) :
  FlexibleStorage1D<T,ST>(toCopy), name_("unfs1d")
{
}

template<typename T, typename ST>
/*virtual*/ const std::string& NamedFlexibleStorage1D<T,ST>::name() const
{
  return name_;
}

template<typename T, typename ST>
void NamedFlexibleStorage1D<T,ST>::operator=(const NamedFlexibleStorage1D<T,ST>& toCopy)
{
  FlexibleStorage1D<T,ST>::operator=(toCopy);
}

template<typename T, typename ST>
void NamedFlexibleStorage1D<T,ST>::operator=(const FlexibleStorage1D<T,ST>& toCopy)
{
  FlexibleStorage1D<T,ST>::operator=(toCopy);
}



#endif