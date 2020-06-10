/**************** written by Thomas Schoenemann as a private person, June 2020 ********************/

#ifndef SORTED_MAP_HH
#define SORTED_MAP_HH

#include <vector>
#include "routines.hh"
#include "flexible_storage1D.hh"
#include "unsorted_map.hh"

//so far operator< is fixed

//use this if your keys are small and easily comparable
//so far this map does not offer erasing
template<typename Key, typename Value, typename KVec = std::vector<Key>, typename VVec = std::vector<Value> >
class SortedMap : public MapBase<Key,Value,KVec,VVec> {
public:

  using Base = MapBase<Key,Value,KVec,VVec>;
  using KeyPassType = typename std::conditional<std::is_fundamental<Key>::value || std::is_pointer<Key>::value, const Key, const Key&>::type;

  SortedMap() {};

  SortedMap(const SortedMap<Key, Value>& toCopy) : MapBase<Key,Value,KVec,VVec>(toCopy) {}
  
  SortedMap(SortedMap<Key, Value>&& toTake) : MapBase<Key,Value,KVec,VVec>(toTake) {}
    
  size_t keypos(KeyPassType key) const noexcept;

  void clear() noexcept { Base::key_.clear(); Base::value_.clear(); }
  
  bool contains(KeyPassType key) noexcept;  

  Value& operator[](KeyPassType key) noexcept;  
};

template<typename Key, typename Value, typename KVec, typename VVec>
size_t SortedMap<Key,Value,KVec,VVec>::keypos(KeyPassType key) const noexcept
{
  return Routines::binsearch_insertpos(Base::key_.data(), key, Base::key_.size());
}

template<typename Key, typename Value, typename KVec, typename VVec>
bool SortedMap<Key,Value,KVec,VVec>::contains(KeyPassType key) noexcept
{
  const size_t size = Base::key_.size();
  return (Routines::binsearch_insertpos(Base::key_.data(), key, size) < size);
}

template<typename Key, typename Value, typename KVec, typename VVec>
Value& SortedMap<Key,Value,KVec,VVec>::operator[](KeyPassType key) noexcept
{ 
  const size_t size = Base::key_.size();
  const size_t inspos = Routines::binsearch_insertpos(Base::key_.data(), key, size);
  
  if (inspos >= size) {
    Base::key_.push_back(key);
    Base::value_.push_back(Value());
    return Base::value_.back();
  }
  else {

    if (Base::key_[inspos] != key) {
            
      Base::key_.push_back(Key());
      Base::value_.push_back(Value());
      Routines::upshift_array(Base::key_.data(), inspos, size, 1);
      Routines::upshift_array(Base::value_.data(), inspos, size, 1);
      Base::key_[inspos] = key;
      Base::value_[inspos] = Value();
    }
    return Base::value_[inspos];
  }
}

#endif