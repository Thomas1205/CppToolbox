/************* written by Thomas Schoenemann as a private person, February 2020 ********/

#ifndef SORTEDSET_HH
#define SORTEDSET_HH

#include "stl_util.hh"
#include "routines.hh"

template<typename T>
class SortedSet {
public:

  using PassType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;

  SortedSet() {}

  SortedSet(const SortedSet<T>& toCopy);
  
  SortedSet(SortedSet<T>&& toTake);
  
  SortedSet(const std::initializer_list<T>& init);

  void swap(SortedSet<T>& other)
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

  //for compatibility with the other sets, e.g. use in templates (data are always sorted)
  const std::vector<T>& unsorted_data() const
  {
    return data_;
  }

  const std::vector<T>& sorted_data() const
  {
    return data_;
  }

  bool contains(PassType val) const;

  //returns true if val is new
  bool insert(PassType val);

  //returns true if val is new
  bool move_insert(T&& val);

  void insert_new(PassType val);

  void move_insert_new(T&& val);
  
  void insert_largest(PassType val);

  void move_insert_largest(T&& val);

  //returns true if val was in the tree
  bool erase(PassType val);

  //returns true if out was in the tree
  bool replace(PassType out, PassType in);

  //returns true if out was in the tree
  bool move_replace(PassType out, T&& in);

protected:

  std::vector<T> data_;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const SortedSet<T>& set);

template<typename T>
bool operator==(const SortedSet<T>& set1, const SortedSet<T>& set2);

/********************** implementation ************************/

template<typename T> 
SortedSet<T>::SortedSet(const SortedSet<T>& toCopy)
  : data_(toCopy.data_) {}

template<typename T> 
SortedSet<T>::SortedSet(SortedSet<T>&& toTake)
 : data_(toTake.data_) {}

template<typename T> 
SortedSet<T>::SortedSet(const std::initializer_list<T>& init)
{
  data_.reserve(init.size());
  for (typename std::initializer_list<T>::const_iterator it = init.begin(); it != init.end(); it++)
    insert(*it); 
}

template<typename T>
bool SortedSet<T>::contains(PassType val) const
{
  return (binsearch(data_, val) != MAX_UINT);
}

//returns true if val is new
template<typename T>
bool SortedSet<T>::insert(PassType val)
{
  //std::cerr << "insert" << std::endl;
  const size_t size = data_.size();
  const size_t inspos = binsearch_insertpos(data_, val);
  if (inspos >= size) {
    data_.push_back(val);
    return true;
  }

  if (data_[inspos] == val)
    return false;

  data_.push_back(T());

  Routines::upshift_array(data_.data(), inspos, size, 1);
  //for (uint k = size; k > inspos; k--)
  //  data_[k] = data_[k-1];

  data_[inspos] = val;
  return true;
}

//returns true if val is new
template<typename T>
bool SortedSet<T>::move_insert(T&& val)
{
  //std::cerr << "insert" << std::endl;
  const size_t size = data_.size();
  const size_t inspos = binsearch_insertpos(data_, val);
  if (inspos >= size) {
    data_.push_back(val);
    return true;
  }

  if (data_[inspos] == val)
    return false;

  data_.push_back(T());

  Routines::upshift_array(data_.data(), inspos, size, 1);
  //for (uint k = size; k > inspos; k--)
  //  data_[k] = data_[k-1];

  data_[inspos] = val;
  return true;
}

//returns true if val is new
template<typename T>
void SortedSet<T>::insert_new(PassType val)
{
  //std::cerr << "insert" << std::endl;
  const size_t size = data_.size();
  const size_t inspos = binsearch_insertpos(data_, val);
  assert(inspos >= size || data_[inspos] != val);
  if (inspos >= size) {
    data_.push_back(val);
  }

  data_.push_back(T());

  Routines::upshift_array(data_.data(), inspos, size, 1);
  //for (uint k = size; k > inspos; k--)
  //  data_[k] = data_[k-1];

  data_[inspos] = val;
}

//returns true if val is new
template<typename T>
void SortedSet<T>::move_insert_new(T&& val)
{
  //std::cerr << "insert" << std::endl;
  const size_t size = data_.size();
  const size_t inspos = binsearch_insertpos(data_, val);
  assert(inspos >= size || data_[inspos] != val);
  if (inspos >= size) {
    data_.push_back(val);
  }

  data_.push_back(T());

  Routines::upshift_array(data_.data(), inspos, size, 1);
  //for (uint k = size; k > inspos; k--)
  //  data_[k] = data_[k-1];

  data_[inspos] = val;
}

template<typename T>
void SortedSet<T>::insert_largest(PassType val)
{
  assert(data_.size() == 0 || data_.back() < val);
  data_.push_back(val);
}

template<typename T>
void SortedSet<T>::move_insert_largest(T&& val)
{
  assert(data_.size() == 0 || data_.back() < val);
  data_.push_back(val);
}

//returns true if val was in the tree
template<typename T>
bool SortedSet<T>::erase(PassType val)
{
  //std::cerr << "erase " << val << " from " << data_ << std::endl;
  const size_t pos = binsearch(data_, val);
  if (pos == MAX_UINT)
    return false;

  const size_t size = data_.size();
  Routines::downshift_array(data_.data(), pos, 1, size);
  //for (uint k = pos; k < size-1; k++)
  //  data_[k] = data_[k+1];

  data_.resize(size-1);
  return true;
}

//returns true if out was in the tree
template<typename T>
bool SortedSet<T>::replace(PassType out, PassType in)
{
  assert(!contains(in));

#if 0
  bool b = erase(out);
  insert(in);
  return b;
#else
  const size_t size = data_.size();
  const size_t pos = binsearch(data_, out);
  if (pos < size) {

    if (pos > 0 && in < data_[pos-1]) {
      size_t npos = pos-1;
      while (npos > 0 && in < data_[npos-1])
        npos--;

      for (size_t k = pos; k > npos; k--)
        data_[k] = data_[k-1];
      data_[npos] = in;
    }
    else if (pos+1 < size && data_[pos+1] < in) {
      size_t npos = pos+1;
      while (npos+1 < size && data_[npos+1] < in)
        npos++;

      for (size_t k = pos; k < npos; k++)
        data_[k] = data_[k+1];
      data_[npos] = in;
    }
    else
      data_[pos] = in;

    return true;
  }
  else {

    insert(in);
    return false;
  }
#endif
}

//returns true if out was in the tree
template<typename T>
bool SortedSet<T>::move_replace(PassType out, T&& in)
{
  assert(!contains(in));

#if 0
  bool b = erase(out);
  insert(in);
  return b;
#else
  const size_t size = data_.size();
  const size_t pos = binsearch(data_, out);
  if (pos < size) {

    if (pos > 0 && in < data_[pos-1]) {
      size_t npos = pos-1;
      while (npos > 0 && in < data_[npos-1])
        npos--;

      for (size_t k = pos; k > npos; k--)
        data_[k] = data_[k-1];
      data_[npos] = in;
    }
    else if (pos+1 < size && data_[pos+1] < in) {
      size_t npos = pos+1;
      while (npos+1 < size && data_[npos+1] < in)
        npos++;

      for (size_t k = pos; k < npos; k++)
        data_[k] = data_[k+1];
      data_[npos] = in;
    }
    else
      data_[pos] = in;

    return true;
  }
  else {

    insert(in);
    return false;
  }
#endif
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const SortedSet<T>& set)
{
  const std::vector<T>& data = set.sorted_data();
  const size_t size = data.size();

  os << "{ ";
  for (size_t i = 0; i < size; i++) {
    if (i != 0)
      os << ", ";
    os << data[i];
  }

  os << " }";
  return os;
}

template<typename T>
bool operator==(const SortedSet<T>& set1, const SortedSet<T>& set2)
{
  return (set1.sorted_data() == set2.sorted_data());
}

#endif