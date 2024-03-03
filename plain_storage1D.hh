/************ written by Thomas Schoenemann, January 2024 *******************/
/*********** Storage1D with no virtual methods (like name), so no pointer to vtable needs to be stored.
             This can be handy when you have a large number of storages ********************/
			 
#ifndef PLAIN_STORAGE1D_HH
#define PLAIN_STORAGE1D_HH

#include "makros.hh"

template<typename T, typename ST=uchar> //if you want to avoid the virtual pointer, most likely your length fits into a uchar
class PlainStorage1D {
public:

  explicit PlainStorage1D() noexcept {}
  
  explicit PlainStorage1D(ST size) : size_(size) { //new throws exceptions
	data_ = new T[size];
  }
  
  explicit PlainStorage1D(ST size, T default_value) : size_(size) {  //new throws exceptions
	data_ = new T[size];
    std::fill_n(data_,size_,default_value);
  }
  
  PlainStorage1D(const std::initializer_list<T>& init) //new throws exceptions
  {
	size_ = init.size();
	data_ = new T[size_];
    std::copy(init.begin(),init.end(),data_);
  }

  PlainStorage1D(const PlainStorage1D<T,ST>& toCopy) { //new throws exceptions
    if (toCopy.data_ == 0) {
	  data_ = 0;
	  size_ = 0;
	  assert(toCopy.size_ == 0);
	}
	else {
  	  size_ = toCopy.size_;
	  data_ = new T[size_];
      Makros::unified_assign(data_, toCopy.direct_access(), size_);
	}
  }
  
  PlainStorage1D(PlainStorage1D<T,ST>&& toTake) noexcept {
	size_ = toTake.size_;
	data_ = toTake.data_;
	toTake.data_ = 0;
  }

  ~PlainStorage1D() {
	delete[] data_;
  }


  PlainStorage1D<T,ST>& operator=(const PlainStorage1D<T,ST>& toCopy)
  {
	delete[] data_;
	if (toCopy.data_ == 0)
	{
	  data_ = 0;
	  size_ = 0;
	  assert(toCopy.size_ == 0);
	}
	else {
  	  size_ = toCopy.size_;
	  data_ = new T[size_];
      Makros::unified_assign(data_, toCopy.direct_access(), size_);	  
	}
	return *this;
  }
  
  PlainStorage1D<T,ST>& operator=(PlainStorage1D<T,ST>&& toTake) noexcept
  {
    delete[] data_;
    data_ = toTake.data_;
    size_ = toTake.size_;
    toTake.data_ = 0;

    return *this;
  }  

  inline const T& operator[](ST i) const noexcept
  {
#ifdef SAFE_MODE
    if (i >= size_) {

      INTERNAL_ERROR << "    invalid const access on element " << ((size_t) i)
                   << " for PlainStorage1D of type "
                   << Makros::Typename<T>()
                   << " with " << ((size_t) size_) << " elements. exiting." << std::endl;

      print_trace();
      exit(1);
    }
#endif
    return data_[i];  
  }

  inline T& operator[](ST i) noexcept
  {
#ifdef SAFE_MODE
    if (i >= size_) {

      INTERNAL_ERROR << "    invalid access on element " << ((size_t) i)
                   << " for PlainStorage1D of type "
                   << Makros::Typename<T>()
                   << " with " << ((size_t) size_) << " elements. exiting." << std::endl;

      print_trace();
      exit(1);
    }
#endif	 
    return data_[i]; 
  }

  void resize(ST new_size)
  {
    if (data_ == 0) {
      data_ = new T[new_size];
    }
    else if (size_ != new_size) {

      T* new_data = new T[new_size];

      const ST size = std::min(size_,new_size);

      Makros::unified_move_assign(new_data, data_, size);

      delete[] data_;
	  data_ = new_data;
    }

    size_ = new_size;  
  }
  
  void resize(ST new_size, T fill_value)
  {
    if (data_ == 0) {
      data_ = new T[new_size];
	  
	  std::fill(data_, data_+new_size, fill_value); //fill and fill_n are of equal speed
    }
    else if (size_ != new_size) {

      T* new_data = new T[new_size];

      const ST size = std::min(size_,new_size);

      Makros::unified_move_assign(new_data, data_, size);

      if (new_size > size_)
        std::fill_n(new_data+size_,new_size-size_,fill_value);

      delete[] data_;
	  data_ = new_data;
    }

    size_ = new_size;  
  }  

  void resize_dirty(ST new_size) { //new throws exceptions

	if (size_ != new_size) {
	  if (data_ != 0)
        delete[] data_;

	  data_ = new T[new_size];
    }
    size_ = new_size;
  }

  inline ST size() const noexcept { return size_; }
  
  inline T* direct_access() noexcept { return data_; }

  inline const T* direct_access() const noexcept { return data_; }
  
  inline const T* end_ptr() const noexcept { return data_ + size_; };

  inline T* end_ptr() noexcept { return data_ + size_; };
  
  T back() const noexcept {
	assert(size_ > 0);
	return data_[size_-1];
  }
  
  inline void set_constant(const T constant) noexcept
  {
	 std::fill_n(data_,size_,constant); 
  }
 

protected:
  T* data_ = 0;	
  ST size_ = 0;
};

template<typename T, typename ST>
bool operator==(const PlainStorage1D<T,ST>& s1, const PlainStorage1D<T,ST>& s2) {
  if (s1.size() != s2.size())
	return false;

  for (uint k=0; k < s1.size(); k++) {
	if (s1.direct_access()[k] != s2.direct_access()[k])
	  return false;
  }
  
  return true;
}

template<typename T, typename ST>
bool operator!=(const PlainStorage1D<T,ST>& s1, const PlainStorage1D<T,ST>& s2) {
  if (s1.size() != s2.size())
	return true;

  for (uint k=0; k < s1.size(); k++) {
	if (s1.direct_access()[k] != s2.direct_access()[k])
	  return true;
  }
  
  return false;
}

template<typename T, typename ST>
bool operator<(const PlainStorage1D<T,ST>& s1, const PlainStorage1D<T,ST>& s2) {

  if (s1.size() != s2.size())
	return (s1.size() < s2.size());

  for (ST k = 0; k < s1.size(); k++) {
	if (s1[k] != s2[k])
	  return (s1[k] < s2[k]);
  }

  return false;
}

template<typename T,typename ST>
std::ostream& operator<<(std::ostream& s, const PlainStorage1D<T,ST>& v)
{
  s << "[ ";
  for (int i=0; i < ((int) v.size()) - 1; i++)
    s << v[i] << ",";
  if (v.size() > 0)
    s << v[v.size()-1];
  s << " ]";

  return s;
}


#endif