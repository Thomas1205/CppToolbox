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
inline typename std::vector<T>::const_iterator set_find(const std::vector<T>& vec, const T val);

template<>
inline std::vector<uint>::const_iterator set_find(const std::vector<uint>& vec, const uint val);

template<typename T>
inline bool set_contains(const std::vector<T>& vec, const T val);

template<>
inline bool set_contains(const std::vector<uint>& vec, const uint val);

template<typename T>
class UnsortedSet {
public:	
  
  using PassType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;
  
  UnsortedSet() {}

  UnsortedSet(const UnsortedSet<T>& toCopy);

  UnsortedSet(UnsortedSet<T>&& toTake);

  UnsortedSet(const std::initializer_list<T>& init);

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

  //const qualifier doesn't make sense here - not returning a reference
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

  bool contains(PassType val) const;

  //returns true if val is new
  bool insert(PassType val);

  //returns true if val is new
  bool move_insert(T&& val);

  void insert_new(PassType val);

  void move_insert_new(T&& val);

  //for compatibility with the other sets (use in templates etc.)
  inline void insert_largest(PassType val)
  {
    insert_new(val);
  }

  //for compatibility with the other sets (use in templates etc.)
  inline void move_insert_largest(T&& val)
  {
    move_insert_new(val);
  }

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
std::ostream& operator<<(std::ostream& os, const UnsortedSet<T>& set);

/********************/

template<typename T>
class UnsortedSetExploitSort {
public:	

  using PassType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;

  UnsortedSetExploitSort() {}

  UnsortedSetExploitSort(const UnsortedSetExploitSort<T>& toCopy);

  UnsortedSetExploitSort(UnsortedSetExploitSort<T>&& toTake);

  UnsortedSetExploitSort(const std::initializer_list<T>& init);

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
      USET_SORT_ALG(data_.data(), data_.size());
      //std::cerr << "result: " << result << std::endl;
      assert(is_unique_sorted(data_.data(), data_.size()));
      is_sorted_ = true;
    }
    
    return data_;
  }

  void get_sorted_data(Storage1D<T>& target) const 
  {
    assign(target, data_);
    USET_SORT_ALG(target.direct_access(), target.size());    
  }

  bool contains(PassType val) const;

  //returns true if val is new
  bool insert(PassType val);

  //returns true if val is new
  bool move_insert(T&& val);

  void insert_new(PassType val);

  void move_insert_new(T&& val);

  //for compatibility with the other sets (use in templates etc.)
  inline void insert_largest(PassType val);

  //for compatibility with the other sets (use in templates etc.)
  inline void move_insert_largest(T&& val);

  //returns true if val was in the tree
  bool erase(PassType val);

  //returns true if out was in the tree
  bool replace(PassType out, PassType in);

  //returns true if out was in the tree
  bool move_replace(PassType out, T&& in);

protected:

  std::vector<T> data_;
  bool is_sorted_ = true;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const UnsortedSetExploitSort<T>& set);

/********************** implementation of helpers ************************/

//find in a sequence without duplicates
template<typename T>
inline typename std::vector<T>::const_iterator set_find(const std::vector<T>& vec, const T val) 
{
  if (std::is_trivially_copyable<T>::value) {
    const uint pos = Routines::find_unique(vec.data(), val, vec.size());
    if (pos >= vec.size())
      return vec.end();
    else
      return vec.begin() + pos;
  }
  else
    return std::find(vec.begin(), vec.end(), val);
}

template<>
inline std::vector<uint>::const_iterator set_find(const std::vector<uint>& vec, const uint val) 
{
  const uint pos = Routines::find_unique_uint(vec.data(), val, vec.size());
  if (pos >= vec.size())
    return vec.end();
  else
    return (vec.begin() + pos);
}

template<typename T>
inline bool set_contains(const std::vector<T>& vec, const T val) 
{
  if (std::is_trivially_copyable<T>::value) 
    return Routines::contains(vec.data(), val);
  else
    return (std::find(vec.begin(), vec.end(), val) != vec.end());
}

template<>
inline bool set_contains(const std::vector<uint>& vec, const uint val) 
{
  return Routines::contains_uint(vec.data(), val, vec.size());
}

/********************** implementation of UnsortedSet ************************/

template<typename T>
UnsortedSet<T>::UnsortedSet(const UnsortedSet<T>& toCopy) 
{
  data_ = toCopy.data_;
}

template<typename T>
UnsortedSet<T>::UnsortedSet(UnsortedSet<T>&& toCopy) 
{
  data_ = toCopy.data_;
}

template<typename T>
UnsortedSet<T>::UnsortedSet(const std::initializer_list<T>& init)
{
  data_.reserve(init.size());
  for (typename std::initializer_list<T>::const_iterator it = init.begin(); it != init.end(); it++)
    insert(*it);
}

template<typename T>
bool UnsortedSet<T>::contains(PassType val) const
{
  //return (set_find(data_, val) != data_.end());
  return set_contains(data_, val);
}

//returns true if val is new
template<typename T>
bool UnsortedSet<T>::insert(PassType val)
{
  //if (set_find(data_, val) != data_.end())
  //  return false;
  if (set_contains(data_, val))
    return false;

  data_.push_back(val);
  return true;
}

//returns true if val is new
template<typename T>
bool UnsortedSet<T>::move_insert(T&& val)
{
  //if (set_find(data_, val) != data_.end())
  //  return false;
  if (set_contains(data_, val))
    return false;

  data_.push_back(val);
  return true;
}

template<typename T>
void UnsortedSet<T>::insert_new(PassType val)
{
  assert(!contains(val));
  data_.push_back(val);
}

template<typename T>
void UnsortedSet<T>::move_insert_new(T&& val)
{
  assert(!contains(val));
  data_.push_back(val);
}

//returns true if val was in the tree
template<typename T>
bool UnsortedSet<T>::erase(PassType val)
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
bool UnsortedSet<T>::replace(PassType out, PassType in)
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

//returns true if out was in the tree
template<typename T>
bool UnsortedSet<T>::move_replace(PassType out, T&& in)
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

template<typename T>
std::ostream& operator<<(std::ostream& os, const UnsortedSet<T>& set)
{
  const std::vector<T> data = set.sorted_data();
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

/********************** implementation of UnsortedSetExploitSort ************************/

template<typename T>
UnsortedSetExploitSort<T>::UnsortedSetExploitSort(const UnsortedSetExploitSort<T>& toCopy)
{
  data_ = toCopy.data_;
  is_sorted_ = toCopy.is_sorted_;
}

template<typename T>
UnsortedSetExploitSort<T>::UnsortedSetExploitSort(UnsortedSetExploitSort<T>&& toTake)
{
  data_ = toTake.data_;
  is_sorted_ = toTake.is_sorted_;
}

template<typename T>
UnsortedSetExploitSort<T>::UnsortedSetExploitSort(const std::initializer_list<T>& init)
{
  data_.reserve(init.size());
  for (typename std::initializer_list<T>::const_iterator it = init.begin(); it != init.end(); it++)
    insert(*it);  
}

template<typename T>
bool UnsortedSetExploitSort<T>::contains(PassType val) const
{
  if (!is_sorted_) {
    //return (set_find(data_, val) != data_.end());
    return set_contains(data_, val);
  }
  else 
    return (Routines::binsearch(data_.data(), val, data_.size()) != MAX_UINT);
}

//returns true if val is new
template<typename T>
bool UnsortedSetExploitSort<T>::insert(PassType val)
{
  const size_t size = data_.size();
  bool is_new = false;
  if (!is_sorted_) {
    //is_new = (set_find(data_, val) != data_.end());
    is_new = !set_contains(data_, val);
  }
  else 
    is_new = (Routines::binsearch(data_.data(), val, size) != MAX_UINT);
  
  if (!is_new) {
    if (is_sorted_ && size > 0 && val < data_.back())
      is_sorted_ = false;
    data_.push_back(val);
  }
  return is_new;
}

//returns true if val is new
template<typename T>
bool UnsortedSetExploitSort<T>::move_insert(T&& val)
{
  const size_t size = data_.size();
  bool is_new = false;
  if (!is_sorted_) {
    //is_new = (set_find(data_, val) != data_.end());
    is_new = !set_contains(data_, val);
  }
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
void UnsortedSetExploitSort<T>::insert_new(PassType val)
{
  assert(!contains(val));
  if (is_sorted_) {
    if (data_.size() > 0 && val < data_.back())
      is_sorted_ = false;
  }
  data_.push_back(val);  
}

template<typename T>
void UnsortedSetExploitSort<T>::move_insert_new(T&& val)
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
inline void UnsortedSetExploitSort<T>::insert_largest(PassType val)
{
  assert(!contains(val));
  data_.push_back(val);  
}

//for compatibility with the other sets (use in templates etc.)
template<typename T>
inline void UnsortedSetExploitSort<T>::move_insert_largest(T&& val)
{
  assert(!contains(val));
  data_.push_back(val);  
}

//returns true if val was in the tree
template<typename T>
bool UnsortedSetExploitSort<T>::erase(PassType val)
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
bool UnsortedSetExploitSort<T>::replace(PassType out, PassType in)
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

//returns true if out was in the tree
template<typename T>
bool UnsortedSetExploitSort<T>::move_replace(PassType out, T&& in)
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

template<typename T>
std::ostream& operator<<(std::ostream& os, const UnsortedSetExploitSort<T>& set)
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

#endif