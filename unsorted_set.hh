/************** written by Thomas Schoenemann as a private person, February 2020 ***********/

#ifndef UNSORTED_SET_HH
#define UNSORTED_SET_HH

#include "sorting.hh"

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
  return (std::find(data_.begin(), data_.end(), val) != data_.end());
}

//returns true if val is new
template<typename T>
bool UnsortedSet<T>::insert(T val)
{
  if (std::find(data_.begin(), data_.end(), val) != data_.end())
    return false;

  data_.push_back(val);
  return true;
}

template<typename T>
void UnsortedSet<T>::insert_new(T val)
{
  assert(std::find(data_.begin(), data_.end(), val) == data_.end());
  data_.push_back(val);
}

//returns true if val was in the tree
template<typename T>
bool UnsortedSet<T>::erase(T val)
{
  typename std::vector<T>::const_iterator it = std::find(data_.begin(), data_.end(), val);
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
  assert(std::find(data_.begin(), data_.end(), in) == data_.end());
	
  typename std::vector<T>::iterator it = std::find(data_.begin(), data_.end(), out);
  if (it == data_.end()) {
    data_.push_back(in);
    return false;
  }
  
  *it = in;
  return true;
}

#endif