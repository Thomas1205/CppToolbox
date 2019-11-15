/*-*-c++-*-*/
/*** written by Thomas Schoenemann as a private person, November 2019 ***/

#ifndef STORGE_BASE_HH
#define STORGE_BASE_HH

#include "makros.hh"

template<typename T, typename ST=size_t>
class StorageBase {
public:

  StorageBase();

  StorageBase(ST size);

  StorageBase(ST size, const T default_value);

	~StorageBase();

  virtual const std::string& name() const;

  inline ST size() const;

  inline T* direct_access();

  inline const T* direct_access() const;

  inline T& direct_access(ST i);

  inline T direct_access(ST i) const;

  inline void set_constant(const T constant);

protected:

  T* data_; //pointers returned by new are guaranteed to have an address that is divisible by 16
  ST size_;
	static const std::string stor_base_name_;
};

/*************************** implementation ****************************/

template<typename T, typename ST>
StorageBase<T,ST>::StorageBase() : data_(0), size_(0)
{
}

template<typename T, typename ST>
StorageBase<T,ST>::StorageBase(ST size) : size_(size)
{
	data_ = new T[size];
}

template<typename T, typename ST>
StorageBase<T,ST>::StorageBase(ST size, const T default_value) : size_(size)
{
	data_ = new T[size];
	std::fill_n(data_,size_,default_value);
}

template<typename T, typename ST>
StorageBase<T,ST>::~StorageBase() 
{
	delete[] data_;
}	

template<typename T, typename ST>
/*virtual*/ const std::string& StorageBase<T,ST>::name() const 
{
	return stor_base_name_;
}

template<typename T, typename ST>
inline ST StorageBase<T,ST>::size() const
{
	return size_;
}

template<typename T, typename ST>
inline T* StorageBase<T,ST>::direct_access()
{
	return data_;
}

template<typename T, typename ST>
inline const T* StorageBase<T,ST>::direct_access() const
{
	return data_;
}

template<typename T, typename ST>
inline T& StorageBase<T,ST>::direct_access(ST i)
{
	return data_[i];
}

template<typename T, typename ST>
inline T StorageBase<T,ST>::direct_access(ST i) const
{
	return data_[i];
}

template<typename T, typename ST>
inline void StorageBase<T,ST>::set_constant(const T constant)
{
	std::fill_n(data_,size_,constant);
}

template<typename T,typename ST>
/*static*/ const std::string StorageBase<T,ST>::stor_base_name_ = "unnamed base storage";


#endif