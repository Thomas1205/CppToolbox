/**************** written by Thomas Schoenemann as a private person, June 2020 ********************/

#ifndef UNSORTED_MAP_HH
#define UNSORTED_MAP_HH

#include <vector>
#include "routines.hh"
#include "flexible_storage1D.hh"

//so far operator< is fixed

template<typename Key, typename Value, typename KVec = std::vector<Key>, typename VVec = std::vector<Value> >
class MapBase {
public: 

  using KeyPassType = typename std::conditional<std::is_fundamental<Key>::value || std::is_pointer<Key>::value, const Key, const Key&>::type;

  MapBase() {};
  
  MapBase(const MapBase<Key, Value>& toCopy) : key_(toCopy.key_), value_(toCopy.value_) {}
  
  MapBase(MapBase<Key, Value>&& toTake) : key_(std::move(toTake.key_)), value_(std::move(toTake.value_)) {}
 
  size_t size() const noexcept { return key_.size(); }

  void reserve(size_t reservedSize) {
    key_.reserve(reservedSize);
    value_.reserve(reservedSize);
  }
 
  void swap(MapBase<Key, Value>& toSwap) noexcept
  {
    key_.swap_(toSwap.key_);
    value_.swap(toSwap.value_);
  }
 
   //NOTE: keys may be unsorted
  const KVec& key() const noexcept { return key_; }
  
  //NOTE: values correspond to keys, which may be unsorted
  const VVec& value() const noexcept { return value_; }
 
protected:

  KVec key_;
  VVec value_;  
};

//use this if your keys are small and easily comparable
//so far this map does not offer erasing
template<typename Key, typename Value, typename KVec = std::vector<Key>, typename VVec = std::vector<Value> >
class UnsortedMap : public MapBase<Key,Value,KVec,VVec> {
public:

  using Base = MapBase<Key,Value,KVec,VVec>;
  using KeyPassType = typename std::conditional<std::is_fundamental<Key>::value || std::is_pointer<Key>::value, const Key, const Key&>::type;

  UnsortedMap() {};
  
  UnsortedMap(const UnsortedMap<Key, Value>& toCopy) : MapBase<Key,Value,KVec,VVec>(toCopy) {}
  
  UnsortedMap(UnsortedMap<Key, Value>&& toTake) : MapBase<Key,Value,KVec,VVec>(toTake) {}
    
  size_t keypos(KeyPassType key) const noexcept;

  void clear() noexcept { Base::key_.clear(); Base::value_.clear(); }
  
  bool contains(KeyPassType key) noexcept;  

  Value& operator[](KeyPassType key) noexcept; 
};

template<typename Key, typename Value, typename KVec = std::vector<Key>, typename VVec = std::vector<Value> >
class UnsortedMapExploitSort : public MapBase<Key,Value,KVec,VVec> {
public:

  using Base = MapBase<Key,Value,KVec,VVec>;
  using KeyPassType = typename std::conditional<std::is_fundamental<Key>::value || std::is_pointer<Key>::value, const Key, const Key&>::type;
  
  UnsortedMapExploitSort() {};
  
  UnsortedMapExploitSort(const UnsortedMapExploitSort<Key, Value>& toCopy) : MapBase<Key,Value,KVec,VVec>(toCopy) {}
  
  UnsortedMapExploitSort(UnsortedMapExploitSort<Key, Value>&& toTake) : MapBase<Key,Value,KVec,VVec>(toTake) {}
  
  void swap(UnsortedMap<Key, Value>& toSwap) noexcept
  {
    Base::key_.swap_(toSwap.key_);
    Base::value_.swap(toSwap.value_);
    std::swap(is_sorted_,toSwap.is_sorted_);
  }
  
  size_t keypos(KeyPassType key) const noexcept;

  void clear() noexcept { Base::key_.clear(); Base::value_.clear(); is_sorted_ = true; }

  bool contains(KeyPassType key) noexcept;  

  Value& operator[](KeyPassType key) noexcept; 
   
protected:

  bool is_sorted_ = true;
};

/****************************** implementation *******************************/

template<typename Key, typename Value, typename KVec, typename VVec>
size_t UnsortedMap<Key,Value,KVec,VVec>::keypos(KeyPassType key) const noexcept
{
  return Routines::find_unique(Base::key_.data(), key, Base::key_.size());
}

template<typename Key, typename Value, typename KVec, typename VVec>
bool UnsortedMap<Key,Value,KVec,VVec>::contains(KeyPassType key) noexcept 
{
  return Routines::contains(Base::key_.data(), key, Base::key_.size());
}

template<typename Key, typename Value, typename KVec, typename VVec>
Value& UnsortedMap<Key,Value,KVec,VVec>::operator[](KeyPassType key) noexcept
{
  const size_t size = Base::key_.size();
  const size_t pos = Routines::find_unique(Base::key_.data(), key, size);
  if (pos < size)
    return Base::value_[pos];
  
  Base::key_.push_back(key);
  Base::value_.push_back(Value());
  return Base::value_.back();
}

/*********************/

template<typename Key, typename Value, typename KVec, typename VVec>
size_t UnsortedMapExploitSort<Key,Value,KVec,VVec>::keypos(KeyPassType key) const noexcept
{
  if (is_sorted_)
    return Routines::binsearch(Base::key_.data(), key, Base::key_.size());
  else
    return Routines::find_unique(Base::key_.data(), key, Base::key_.size());
}

template<typename Key, typename Value, typename KVec, typename VVec>
bool UnsortedMapExploitSort<Key,Value,KVec,VVec>::contains(KeyPassType key) noexcept 
{
  const size_t size = Base::key_.size();
  if (is_sorted_)
    return (Routines::binsearch(Base::key_.data(), key, size) < size);
  else
    return Routines::contains(Base::key_.data(), key, size);
}

template<typename Key, typename Value, typename KVec, typename VVec>
Value& UnsortedMapExploitSort<Key,Value,KVec,VVec>::operator[](KeyPassType key) noexcept
{  
  const size_t size = Base::key_.size();
  const size_t pos = (is_sorted_) ? Routines::binsearch(Base::key_.data(), key, size) : Routines::find_unique(Base::key_.data(), key, size);
  if (pos < size)
    return Base::value_[pos];

  if (size > 0 && key < Base::key_.back())
    is_sorted_ = false;

  Base::key_.push_back(key);
  Base::value_.push_back(Value());
  return Base::value_.back();
}

#endif