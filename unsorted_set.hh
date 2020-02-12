/************** written by Thomas Schoenemann as a private person, February 2020 ***********/

#ifndef UNSORTED_SET_HH
#define UNSORTED_SET_HH

#include "sorting.hh"
#include "routines.hh"

#ifndef USET_SORT_ALG
#define USET_SORT_ALG bubble_sort
#endif

//find in a sequence without duplicates
template<typename T>
inline typename std::vector<T>::const_iterator set_find(const std::vector<T>& vec, T val) {
  return std::find(vec.begin(), vec.end(), val);
}

template<>
inline std::vector<uint>::const_iterator set_find(const std::vector<uint>& vec, uint val) {
  const uint pos = Routines::find_unique_uint(vec.data(), val, vec.size());
  if (pos >= vec.size())
    return vec.end();
  else
    return (vec.begin() + pos);
}

template<>
inline std::vector<int>::const_iterator set_find(const std::vector<int>& vec, int val) {
  const uint pos = Routines::find_unique_int(vec.data(), val, vec.size());
  if (pos >= vec.size())
    return vec.end();
  else
    return (vec.begin() + pos);
}

template<typename T>
class UnsortedSet {
public:	
  
  UnsortedSet() {}

  UnsortedSet(const UnsortedSet<T>& toCopy);

  void swap(UnsortedSet<T>& other)
  {
    data_.swap(other.data_);
  }

  size_t size() const
  {
    return data_.size();
  }

  size_t capacity() const
  {
    return data_.capacity();
  }

  void reserve(size_t size)
  {
    data_.reserve(size);
  }

  void clear()
  {
    data_.clear();
  }

  const std::vector<T>& unsorted_data() const
  {
    return data_;
  }

  const std::vector<T> sorted_data() const
  {
    std::vector<T> result = data_;
    USET_SORT_ALG(result.data(), result.size());
    assert(is_unique_sorted(result.data(), result.size()));
    
    return result;
  }

  void get_sorted_data(Storage1D<T>& target) const 
  {
    assign(target, data_);
    USET_SORT_ALG(target.direct_access(), target.size());    
  }

  bool contains(T val) const;

  //returns true if val is new
  bool insert(T val);

  void insert_new(T val);

  //for compatibility with the other sets (use in templates etc.)
  inline void insert_largest(const T val)
  {
    insert_new(val);
  }

  //returns true if val was in the tree
  bool erase(T val);

  //returns true if out was in the tree
  bool replace(T out, T in);

protected:

  std::vector<T> data_;
};

/********************/

template<typename T>
class UnsortedSetExploitSort {
public:	

  UnsortedSetExploitSort() {}

  UnsortedSetExploitSort(const UnsortedSetExploitSort<T>& toCopy);

  void swap(UnsortedSetExploitSort<T>& other)
  {
    data_.swap(other.data_);
    std::swap(is_sorted_, other.is_sorted_);
  }

  size_t size() const
  {
    return data_.size();
  }

  size_t capacity() const
  {
    return data_.capacity();
  }

  void reserve(size_t size)
  {
    data_.reserve(size);
  }

  void clear()
  {
    data_.clear();
  }

  const std::vector<T>& unsorted_data() const
  {
    return data_;
  }

  const std::vector<T>& sorted_data()
  {
    if (!is_sorted_) {
      USET_SORT_ALG(data_, data_.size());
      //std::cerr << "result: " << result << std::endl;
      assert(is_unique_sorted(data_, data_.size()));
      is_sorted_ = true;
    }
    
    return data_;
  }

  void get_sorted_data(Storage1D<T>& target) const 
  {
    assign(target, data_);
    USET_SORT_ALG(target.direct_access(), target.size());    
  }

  bool contains(T val) const;

  //returns true if val is new
  bool insert(T val);

  void insert_new(T val);

  //for compatibility with the other sets (use in templates etc.)
  inline void insert_largest(const T val);

  //returns true if val was in the tree
  bool erase(T val);

  //returns true if out was in the tree
  bool replace(T out, T in);

protected:

  std::vector<T> data_;
  bool is_sorted_ = true;
};

/********************** implementation of UnsortedSet ************************/

template<typename T>
UnsortedSet<T>::UnsortedSet(const UnsortedSet<T>& toCopy) {
  data_ = toCopy.data_;
}

template<typename T>
bool UnsortedSet<T>::contains(T val) const
{
  return (set_find(data_, val) != data_.end());
}

//returns true if val is new
template<typename T>
bool UnsortedSet<T>::insert(T val)
{
  if (set_find(data_, val) != data_.end())
    return false;

  data_.push_back(val);
  return true;
}

template<typename T>
void UnsortedSet<T>::insert_new(T val)
{
  assert(set_find(data_, val) == data_.end());
  data_.push_back(val);
}

//returns true if val was in the tree
template<typename T>
bool UnsortedSet<T>::erase(T val)
{
  const typename std::vector<T>::const_iterator it = set_find(data_, val);
  if (it == data_.end())
    return false;

  size_t pos = it - data_.begin();
  data_[pos] = data_.back();
  data_.resize(data_.size()-1);
  return true;
}

//returns true if out was in the tree
template<typename T>
bool UnsortedSet<T>::replace(T out, T in)
{
  assert(!contains(in));
	
  typename std::vector<T>::const_iterator it = set_find(data_, out);
  if (it == data_.end()) {
    data_.push_back(in);
    return false;
  }
  
  data_[it - data_.begin()] = in;
  return true;
}


/********************** implementation of UnsortedSetExploitSort ************************/

template<typename T>
UnsortedSetExploitSort<T>::UnsortedSetExploitSort(const UnsortedSetExploitSort<T>& toCopy)
{
  data_ = toCopy.data_;
  is_sorted_ = toCopy.is_sorted_;
}

template<typename T>
bool UnsortedSetExploitSort<T>::contains(T val) const
{
  if (!is_sorted_)
    return (set_find(data_, val) != data_.end());
  else 
    return (Routines::binsearch(data_.data(), val, data_.size()) != MAX_UINT);
}

//returns true if val is new
template<typename T>
bool UnsortedSetExploitSort<T>::insert(T val)
{
  const size_t size = data_.size();
  bool is_new = false;
  if (!is_sorted_)
    is_new = (set_find(data_, val) != data_.end());
  else 
    is_new = (Routines::binsearch(data_.data(), val, size) != MAX_UINT);
  
  if (!is_new) {
    if (is_sorted_ && size > 0 && val < data_.back())
      is_sorted_ = false;
    data_.push_back(val);
  }
  return is_new;
}

template<typename T>
void UnsortedSetExploitSort<T>::insert_new(T val)
{
  assert(!contains(val));
  if (is_sorted_) {
    if (data_.size() > 0 && val < data_.back())
      is_sorted_ = false;
  }
  data_.push_back(val);  
}

//for compatibility with the other sets (use in templates etc.)
template<typename T>
inline void UnsortedSetExploitSort<T>::insert_largest(const T val)
{
  data_.push_back(val);  
}

//returns true if val was in the tree
template<typename T>
bool UnsortedSetExploitSort<T>::erase(T val)
{
  const size_t size = data_.size();
  size_t pos = 0;
  if (!is_sorted_) {
    const typename std::vector<T>::const_iterator it = set_find(data_, val);
    if (it == data_.end())
      return false;
    pos = it - data_.begin();
  }
  else {
    pos = Routines::binsearch(data_.data(), val, data_.size());
    if (pos == MAX_UINT)
      return false;
  }
  
  if (pos != size - 1) {
    data_[pos] = data_.back();
    is_sorted_ = false;
  }
  data_.resize(size-1);
  return true;  
}

//returns true if out was in the tree
template<typename T>
bool UnsortedSetExploitSort<T>::replace(T out, T in)
{
  assert(!contains(in));

  const size_t size = data_.size();
  size_t pos = 0;
  if (!is_sorted_) {
    const typename std::vector<T>::const_iterator it = set_find(data_, out);
    if (it == data_.end())
      pos = MAX_UINT;
    else
      pos = it - data_.begin();
  }
  else {
    pos = Routines::binsearch(data_.data(), out, size);
  }
    
  if (pos == MAX_UINT) {
    if (is_sorted_ && size > 0 && in < data_.back())
      is_sorted_ = false;
    data_.push_back(in);
    return false;
  }
  
  is_sorted_ = false;
  data_[pos] = in;
  return true;
}


#endif