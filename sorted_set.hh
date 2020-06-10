/************* written by Thomas Schoenemann as a private person, February 2020 ********/

#ifndef SORTEDSET_HH
#define SORTEDSET_HH

#include "stl_util.hh"
#include "routines.hh"
#include "unsorted_set.hh"

template<typename T>
class SortedSet : public SetBase<T> {
public:

  using Base = SetBase<T>;
  using PassType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;

  SortedSet() {}

  SortedSet(const SortedSet<T>& toCopy) : SetBase<T>(toCopy) {}
  
  SortedSet(SortedSet<T>&& toTake) : SetBase<T>(toTake) {}
  
  SortedSet(const std::initializer_list<T>& init);

  //for compatibility with the other sets, e.g. use in templates (data are always sorted)
  const std::vector<T>& unsorted_data() const noexcept
  {
    return Base::data_;
  }

  const std::vector<T>& sorted_data() const noexcept
  {
    return Base::data_;
  }

  bool contains(PassType val) const noexcept;

  //returns true if val is new
  bool insert(PassType val) noexcept;

  //returns true if val is new
  bool move_insert(T&& val) noexcept;

  void insert_new(PassType val) noexcept;

  void move_insert_new(T&& val) noexcept;
  
  void insert_largest(PassType val) noexcept;

  void move_insert_largest(T&& val) noexcept;

  //returns true if val was in the tree
  bool erase(PassType val) noexcept;

  //returns true if out was in the tree
  bool replace(PassType out, PassType in) noexcept;

  //returns true if out was in the tree
  bool move_replace(PassType out, T&& in) noexcept;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const SortedSet<T>& set);

template<typename T>
bool operator==(const SortedSet<T>& set1, const SortedSet<T>& set2) noexcept;

/********************** implementation ************************/

template<typename T> 
SortedSet<T>::SortedSet(const std::initializer_list<T>& init)
{
  Base::data_.reserve(init.size());
  for (typename std::initializer_list<T>::const_iterator it = init.begin(); it != init.end(); it++)
    insert(*it); 
}

template<typename T>
bool SortedSet<T>::contains(PassType val) const noexcept
{
  return (binsearch(Base::data_, val) < Base::data_.size());
}

//returns true if val is new
template<typename T>
bool SortedSet<T>::insert(PassType val) noexcept
{
  //std::cerr << "insert" << std::endl;
  const size_t size = Base::data_.size();
  const size_t inspos = binsearch_insertpos(Base::data_, val);
  if (inspos >= size) {
    Base::data_.push_back(val);
    return true;
  }

  if (Base::data_[inspos] == val)
    return false;

  Base::data_.push_back(T());

  Routines::upshift_array(Base::data_.data(), inspos, size, 1);
  //for (uint k = size; k > inspos; k--)
  //  data_[k] = data_[k-1];

  Base::data_[inspos] = val;
  return true;
}

//returns true if val is new
template<typename T>
bool SortedSet<T>::move_insert(T&& val) noexcept
{
  //std::cerr << "insert" << std::endl;
  const size_t size = Base::data_.size();
  const size_t inspos = binsearch_insertpos(Base::data_, val);
  if (inspos >= size) {
    Base::data_.push_back(val);
    return true;
  }

  if (Base::data_[inspos] == val)
    return false;

  Base::data_.push_back(T());

  Routines::upshift_array(Base::data_.data(), inspos, size, 1);
  //for (uint k = size; k > inspos; k--)
  //  Base::data_[k] = Base::data_[k-1];

  Base::data_[inspos] = val;
  return true;
}

//returns true if val is new
template<typename T>
void SortedSet<T>::insert_new(PassType val) noexcept
{
  //std::cerr << "insert" << std::endl;
  const size_t size = Base::data_.size();
  const size_t inspos = binsearch_insertpos(Base::data_, val);
  assert(inspos >= size || Base::data_[inspos] != val);
  if (inspos >= size) {
    Base::data_.push_back(val);
  }

  Base::data_.push_back(T());

  Routines::upshift_array(Base::data_.data(), inspos, size, 1);
  //for (uint k = size; k > inspos; k--)
  //  Base::data_[k] = Base::data_[k-1];

  Base::data_[inspos] = val;
}

//returns true if val is new
template<typename T>
void SortedSet<T>::move_insert_new(T&& val) noexcept
{
  //std::cerr << "insert" << std::endl;
  const size_t size = Base::data_.size();
  const size_t inspos = binsearch_insertpos(Base::data_, val);
  assert(inspos >= size || Base::data_[inspos] != val);
  if (inspos >= size) {
    Base::data_.push_back(val);
  }

  Base::data_.push_back(T());

  Routines::upshift_array(Base::data_.data(), inspos, size, 1);
  //for (uint k = size; k > inspos; k--)
  //  Base::data_[k] = Base::data_[k-1];

  Base::data_[inspos] = val;
}

template<typename T>
void SortedSet<T>::insert_largest(PassType val) noexcept
{
  assert(Base::data_.size() == 0 || Base::data_.back() < val);
  Base::data_.push_back(val);
}

template<typename T>
void SortedSet<T>::move_insert_largest(T&& val) noexcept
{
  assert(Base::data_.size() == 0 || Base::data_.back() < val);
  Base::data_.push_back(val);
}

//returns true if val was in the tree
template<typename T>
bool SortedSet<T>::erase(PassType val) noexcept
{
  //std::cerr << "erase " << val << " from " << data_ << std::endl;
  const size_t size = Base::data_.size();
  const size_t pos = binsearch(Base::data_, val);
  if (pos >= size)
    return false;

  Routines::downshift_array(Base::data_.data(), pos, 1, size);
  //for (uint k = pos; k < size-1; k++)
  //  data_[k] = data_[k+1];

  Base::data_.resize(size-1);
  return true;
}

//returns true if out was in the tree
template<typename T>
bool SortedSet<T>::replace(PassType out, PassType in) noexcept
{
  assert(!contains(in));

  //std::cerr << "replace " << out << " by " << in << std::endl;

#if 0
  bool b = erase(out);
  insert(in);
  return b;
#else
  const size_t size = Base::data_.size();
  const size_t pos = binsearch(Base::data_, out);
  if (pos < size) {

    if (pos > 0 && in < Base::data_[pos-1]) {
      size_t npos = pos-1;
      while (npos > 0 && in < Base::data_[npos-1])
        npos--;

      //std::cerr << "A. npos: " << npos << std::endl;

      Routines::upshift_array(Base::data_.data(), npos, pos, 1);
      //for (size_t k = pos; k > npos; k--)
      //  Base::data_[k] = Base::data_[k-1];
      Base::data_[npos] = in;
    }
    else if (pos+1 < size && Base::data_[pos+1] < in) {
      size_t npos = pos+1;
      while (npos+1 < size && Base::data_[npos+1] < in)
        npos++;

      //std::cerr << "B. npos: " << npos << std::endl;

      Routines::downshift_array(Base::data_.data(), pos, 1, npos+1);      
      //for (size_t k = pos; k < npos; k++)
      //  Base::data_[k] = Base::data_[k+1];
      Base::data_[npos] = in;
    }
    else
      Base::data_[pos] = in;

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
bool SortedSet<T>::move_replace(PassType out, T&& in) noexcept
{
  assert(!contains(in));

#if 0
  bool b = erase(out);
  insert(in);
  return b;
#else
  const size_t size = Base::data_.size();
  const size_t pos = binsearch(Base::data_, out);
  if (pos < size) {

    if (pos > 0 && in < Base::data_[pos-1]) {
      size_t npos = pos-1;
      while (npos > 0 && in < Base::data_[npos-1])
        npos--;

      Routines::upshift_array(Base::data_.data(), npos, pos, 1);
      //for (size_t k = pos; k > npos; k--)
      //  Base::data_[k] = Base::data_[k-1];
      Base::data_[npos] = in;
    }
    else if (pos+1 < size && Base::data_[pos+1] < in) {
      size_t npos = pos+1;
      while (npos+1 < size && Base::data_[npos+1] < in)
        npos++;

      Routines::downshift_array(Base::data_.data(), pos, 1, npos+1);      
      //for (size_t k = pos; k < npos; k++)
      //  Base::data_[k] = Base::data_[k+1];
      Base::data_[npos] = in;
    }
    else
      Base::data_[pos] = in;

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
bool operator==(const SortedSet<T>& set1, const SortedSet<T>& set2) noexcept
{
  return (set1.sorted_data() == set2.sorted_data());
}

#endif