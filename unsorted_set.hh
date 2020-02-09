/************** written by Thomas Schoenemann as a private person, February 2020 ***********/

#ifndef UNSORTED_SET_HH
#define UNSORTED_SET_HH

#include "sorting.hh"

//find in a sequence without duplicates
template<typename T>
inline typename std::vector<T>::const_iterator set_find(const std::vector<T>& vec, T val) {
  return std::find(vec.begin(), vec.end(), val);
}

template<>
inline std::vector<uint>::const_iterator set_find(const std::vector<uint>& vec, uint val) {
  const uint pos = Makros::find_unique_uint(vec.data(), val, vec.size());
  if (pos >= vec.size())
    return vec.end();
  else
    return (vec.begin() + pos);
}

template<>
inline std::vector<int>::const_iterator set_find(const std::vector<int>& vec, int val) {
  const uint pos = Makros::find_unique_int(vec.data(), val, vec.size());
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
    bubble_sort(result.data(), result.size());
    //std::cerr << "result: " << result << std::endl;
    assert(is_unique_sorted(result.data(), result.size()));
    
    return result;
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

/********************** implementation ************************/

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
  assert(set_find(data_, in) == data_.end());
	
  typename std::vector<T>::const_iterator it = set_find(data_, out);
  if (it == data_.end()) {
    data_.push_back(in);
    return false;
  }
  
  data_[it - data_.begin()] = in;
  return true;
}

#endif