/********************** written by Thomas Schoenemann as a private person, June 2020 *************************/

#ifndef HASHMAP_HH
#define HASHMAP_HH

#include "flexible_storage1D.hh"
#include "routines.hh"

template<typename Key, typename Value, typename Hash>
class HashMapBase {
public:

  using KeyPassType = typename std::conditional<std::is_fundamental<Key>::value || std::is_pointer<Key>::value, const Key, const Key&>::type;

  HashMapBase() {}
  
  HashMapBase(const HashMapBase& toCopy) : hash_value_(toCopy.hash_value_), key_(toCopy.key_), value_(toCopy.value_) {}

  HashMapBase(HashMapBase&& toTake) : hash_value_(std::move(toTake.hash_value_)), key_(std::move(toTake.key_)), value_(std::move(toTake.value_)) {}

  void swap(HashMapBase& toSwap) {
    hash_value_.swap(toSwap.hash_value_);
    key_.swap(toSwap.key_);
    value_.swap(toSwap.value_);
  }

  struct iterator;
  struct const_iterator;

  bool empty() const noexcept;  
    
  size_t size() const noexcept;

  size_t bucket_count() const noexcept;

  iterator begin() noexcept;
  
  iterator end() noexcept;

  const_iterator cbegin() const noexcept;
  
  const_iterator cend() const noexcept;

  void clear() {
    hash_value_.clear(true);
    key_.clear(true);
    value_.clear(true);
  }

  struct iterator {
    
    iterator(HashMapBase<Key,Value,Hash>& map) : map_(map) {}

    iterator(HashMapBase<Key,Value,Hash>& map, size_t pos, size_t sub_pos) : map_(map), pos_(pos), sub_pos_(sub_pos) {}
    
    inline iterator& operator++() noexcept {      
      
      const FlexibleStorage1D<FlexibleStorage1D<Key> >& map_key = map_.key_;
      
      assert(pos_ < map_key.size());
      if (sub_pos_ + 1 < map_key[pos_].size())
        sub_pos_++;
      else {
        pos_++;
        sub_pos_ = 0;
        if (pos_ < map_key.size()) 
          assert(map_key[pos_].size() > 0);
      }
      
      return *this;
    }
    
    iterator operator++(int) noexcept {
      iterator res = *this;
      operator++();
      return res;
    }

    bool operator!=(const HashMapBase<Key,Value,Hash>::iterator& i2) const noexcept
    {
      return (pos_ != i2.pos_ || sub_pos_ != i2.sub_pos_);      
    }
    
    //no ideas for this yet
    //std::pair<KeyPassType, Value&> operator*();
    //std::pair<KeyPassType, Value&> operator->();
      
    const Key& key() const noexcept {
      return map_.key_[pos_][sub_pos_];
    }
      
    Value& value() {
      return map_.value_[pos_][sub_pos_];
    }
      
    size_t pos_ = 0;
    size_t sub_pos_ = 0;
    HashMapBase<Key,Value,Hash>& map_;
  };

  struct const_iterator {

    const_iterator(const HashMapBase<Key,Value,Hash>& map) : map_(map) {}

    const_iterator(const HashMapBase<Key,Value,Hash>& map, size_t pos, size_t sub_pos) : map_(map), pos_(pos), sub_pos_(sub_pos) {}
        
    bool operator!=(const HashMapBase<Key,Value,Hash>::const_iterator& i2) const noexcept
    {
      return (pos_ != i2.pos_ || sub_pos_ != i2.sub_pos_);      
    }

    bool operator!=(const HashMapBase<Key,Value,Hash>::iterator& i2) const noexcept
    {
      return (pos_ != i2.pos_ || sub_pos_ != i2.sub_pos_);      
    }
          
    inline const_iterator& operator++() noexcept {      
      
      const FlexibleStorage1D<FlexibleStorage1D<Key> >& map_key = map_.key_;
            
      assert(pos_ < map_key.size());
      if (sub_pos_ + 1 < map_key[pos_].size())
        sub_pos_++;
      else {
        pos_++;
        sub_pos_ = 0;
        if (pos_ < map_key.size()) 
          assert(map_key[pos_].size() > 0);
      }
      
      return *this;
    }
    
    const_iterator operator++(int) noexcept {
      const_iterator res = *this;
      operator++();
      return res;
    }
                             
    const Key& key() const noexcept {
      return map_.key_[pos_][sub_pos_];
    }
      
    const Value& value() const noexcept {
      return map_.value_[pos_][sub_pos_];
    }

    size_t pos_ = 0;
    size_t sub_pos_ = 0;
    const HashMapBase<Key,Value,Hash>& map_;    
  };

  inline Value& inner_op(KeyPassType key, size_t pos) {
    
    FlexibleStorage1D<Key>& cur_key_list = key_[pos];
    FlexibleStorage1D<Value>& cur_value_list = value_[pos];
    const size_t sub_pos = Routines::find_unique(cur_key_list.direct_access(), key, cur_key_list.size());
    if (sub_pos < cur_key_list.size())
      return cur_value_list[sub_pos];
    else {    
      cur_key_list.push_back(key);
      return cur_value_list.increase();
    }
  }

protected:

  FlexibleStorage1D<size_t> hash_value_;
  FlexibleStorage1D<FlexibleStorage1D<Key> > key_;
  FlexibleStorage1D<FlexibleStorage1D<Value> > value_;

  static Hash hash_;
};

template<typename TK, typename TV, typename H>
std::ostream& operator<<(std::ostream& os, const HashMapBase<TK,TV,H>& m);


//so far this map does not offer erasing
template<typename Key, typename Value, typename Hash>
class UnsortedHashMap : public HashMapBase<Key, Value, Hash> {
public:

  using KeyPassType = typename std::conditional<std::is_fundamental<Key>::value || std::is_pointer<Key>::value, const Key, const Key&>::type;
  using Base = HashMapBase<Key, Value, Hash>;

  UnsortedHashMap() : HashMapBase<Key,Value,Hash>() {}
  
  UnsortedHashMap(const UnsortedHashMap<Key,Value,Hash>& toCopy) : HashMapBase<Key,Value,Hash>(toCopy) {}

  UnsortedHashMap(UnsortedHashMap<Key,Value,Hash>&& toTake) : HashMapBase<Key,Value,Hash>(toTake) {}
  
  Value& operator[](KeyPassType key); 

  //Value& operator[](Key&& key); 

  //bool try_get(KeyPassType key, Value*& ptr) const;

  bool contains(KeyPassType key);
};

//so far this map does not offer erasing
template<typename Key, typename Value, typename Hash>
class UnsortedHashMapExploitSort : public HashMapBase<Key, Value, Hash> {
public:

  using KeyPassType = typename std::conditional<std::is_fundamental<Key>::value || std::is_pointer<Key>::value, const Key, const Key&>::type;
  using Base = HashMapBase<Key, Value, Hash>;

  UnsortedHashMapExploitSort() : HashMapBase<Key,Value,Hash>() {}
  
  UnsortedHashMapExploitSort(const UnsortedHashMapExploitSort<Key,Value,Hash>& toCopy) : HashMapBase<Key,Value,Hash>(toCopy) {}

  UnsortedHashMapExploitSort(UnsortedHashMapExploitSort<Key,Value,Hash>&& toTake) : HashMapBase<Key,Value,Hash>(toTake) {}
  
  Value& operator[](KeyPassType key); 

  //Value& operator[](Key&& key); 

  //bool try_get(KeyPassType key, Value*& ptr) const;

  bool contains(KeyPassType key);

  void swap(UnsortedHashMapExploitSort& toSwap) 
  {
    Base::swap(toSwap);
    std::swap(is_sorted_, toSwap.is_sorted_);
  }

  void clear() 
  {
    Base::clear();
    is_sorted_ = true;
  }
  
protected: 

  bool is_sorted_ = true;
};

//so far this map does not offer erasing
template<typename Key, typename Value, typename Hash>
class SortedHashMap : public HashMapBase<Key, Value, Hash> {
public:

  using KeyPassType = typename std::conditional<std::is_fundamental<Key>::value || std::is_pointer<Key>::value, const Key, const Key&>::type;
  using Base = HashMapBase<Key, Value, Hash>;

  SortedHashMap() : HashMapBase<Key,Value,Hash>() {}
  
  SortedHashMap(const SortedHashMap<Key,Value,Hash>& toCopy) : HashMapBase<Key,Value,Hash>(toCopy) {}

  SortedHashMap(SortedHashMap<Key,Value,Hash>&& toTake) : HashMapBase<Key,Value,Hash>(toTake) {}

  Value& operator[](KeyPassType key); 

  //Value& operator[](Key&& key); 

  //bool try_get(KeyPassType key, Value*& ptr) const;

  bool contains(KeyPassType key);
};

/******************** implementation **************************/

template<typename TK, typename TV, typename H>
std::ostream& operator<<(std::ostream& os, const HashMapBase<TK,TV,H>& m)
{  
  os << "[ ";
  for (typename HashMapBase<TK,TV,H>::const_iterator it = m.cbegin(); it != m.cend(); ) {
    os << it.key() << "->" << it.value();
    ++it;
    if (it != m.cend())
      os << ", ";
  }
  os << " ]";

  return os;
}

template<typename Key, typename Value, typename Hash>
bool HashMapBase<Key,Value,Hash>::empty() const noexcept
{
  return (key_.size() == 0); 
}

template<typename Key, typename Value, typename Hash>
size_t HashMapBase<Key,Value,Hash>::size() const noexcept
{
  size_t size = 0;
  for (size_t i=0; i < key_.size(); i++)
    size += key_[i].size();

  return size;
}

template<typename Key, typename Value, typename Hash>
size_t HashMapBase<Key,Value,Hash>::bucket_count() const noexcept
{
  return key_.size();
}

template<typename Key, typename Value, typename Hash>
typename HashMapBase<Key,Value,Hash>::iterator HashMapBase<Key,Value,Hash>::begin() noexcept
{
  return iterator(*this);
}
  
template<typename Key, typename Value, typename Hash>
typename HashMapBase<Key,Value,Hash>::iterator HashMapBase<Key,Value,Hash>::end() noexcept
{
  if (key_.size() == 0)
    return iterator(*this);
  else
    return iterator(*this,key_.size(), 0);
}

template<typename Key, typename Value, typename Hash>
typename HashMapBase<Key,Value,Hash>::const_iterator HashMapBase<Key,Value,Hash>::cbegin() const noexcept
{
  //std::cerr << "key: " << key_ << std::endl; 
  return const_iterator(*this);  
}
  
template<typename Key, typename Value, typename Hash>
typename HashMapBase<Key,Value,Hash>::const_iterator HashMapBase<Key,Value,Hash>::cend() const noexcept
{
  if (key_.size() == 0)
    return const_iterator(*this);
  else
    return const_iterator(*this,key_.size(), 0);  
}

/*************/

template<typename Key, typename Value, typename Hash>
Value& UnsortedHashMap<Key,Value,Hash>::operator[](KeyPassType key)
{
  const size_t size = Base::hash_value_.size();
  const size_t cur_hash = Base::hash_(key);
  const size_t pos = Routines::find_unique(Base::hash_value_.direct_access(), cur_hash, size); 
  
  if (pos >= size) {
    Base::hash_value_.append(cur_hash);
    Base::key_.increase().append(key);
    return Base::value_.increase().increase();
  }

  return Base::inner_op(key, pos);
}

template<typename Key, typename Value, typename Hash>
bool UnsortedHashMap<Key,Value,Hash>::contains(KeyPassType key) 
{
  const size_t size = Base::hash_value_.size();
  const size_t cur_hash = Base::hash_(key);
  const size_t pos = Routines::find_unique(Base::hash_value_.direct_access(), cur_hash, size); 
  
  if (pos >= size) 
    return false;
  
  const FlexibleStorage1D<Key>& cur_key_list = Base::key_[pos];
  const size_t sub_pos = Routines::find_unique(cur_key_list.direct_access(), key, cur_key_list.size());
  return (sub_pos < cur_key_list.size());
}

/*********************/

template<typename Key, typename Value, typename Hash>
Value& UnsortedHashMapExploitSort<Key,Value,Hash>::operator[](KeyPassType key)
{
  const size_t size = Base::hash_value_.size();
  const size_t cur_hash = Base::hash_(key);
  const size_t pos = (is_sorted_) ?
    Routines::binsearch(Base::hash_value_.direct_access(), cur_hash, size) :   
    Routines::find_unique(Base::hash_value_.direct_access(), cur_hash, size); 

  if (pos >= size) {
    if (size > 0 && cur_hash < Base::hash_value_.back())
      is_sorted_ = false;
    
    Base::hash_value_.append(cur_hash);
    Base::key_.increase().append(key);
    return Base::value_.increase().increase();
  }

  return Base::inner_op(key, pos);
}

template<typename Key, typename Value, typename Hash>
bool UnsortedHashMapExploitSort<Key,Value,Hash>::contains(KeyPassType key) 
{
  const size_t size = Base::hash_value_.size();
  const size_t cur_hash = Base::hash_(key);
  const size_t pos = (is_sorted_) ?
    Routines::binsearch(Base::hash_value_.direct_access(), cur_hash, size) :   
    Routines::find_unique(Base::hash_value_.direct_access(), cur_hash, size); 
  
  if (pos >= size) 
    return false;
  
  const FlexibleStorage1D<Key>& cur_key_list = Base::key_[pos];
  const size_t sub_pos = Routines::find_unique(cur_key_list.direct_access(), key, cur_key_list.size());
  return (sub_pos < cur_key_list.size());
}

/*********************/

template<typename Key, typename Value, typename Hash>
Value& SortedHashMap<Key,Value,Hash>::operator[](KeyPassType key)
{  
  const size_t size = Base::hash_value_.size();
  const size_t cur_hash = Base::hash_(key);
  const size_t pos = Routines::binsearch_insertpos(Base::hash_value_.direct_access(), cur_hash, size);

  assert(Base::key_.size() == size);
  assert(Base::value_.size() == size);

  if (pos >= size) {
    Base::hash_value_.append(cur_hash);
    Base::key_.increase().append(key);
    return Base::value_.increase().increase();
  }
  else if (Base::hash_value_[pos] != cur_hash) {
    
    Base::hash_value_.increase();
    Base::key_.increase();
    Base::value_.increase();
    
    Routines::upshift_array(Base::hash_value_.data(), pos, size, 1);
    Routines::upshift_array(Base::key_.data(), pos, size, 1);
    Routines::upshift_array(Base::value_.data(), pos, size, 1);
        
    Base::hash_value_[pos] = cur_hash;
    Base::key_[pos] = FlexibleStorage1D<Key>();
    Base::value_[pos] = FlexibleStorage1D<Value>();
    
    Base::key_[pos].append(key);
    return Base::value_[pos].increase();
  }
  
  return Base::inner_op(key, pos);
}

template<typename Key, typename Value, typename Hash>
bool SortedHashMap<Key,Value,Hash>::contains(KeyPassType key) 
{
  const size_t size = Base::hash_value_.size();
  const size_t cur_hash = Base::hash_(key);
  const size_t pos = Routines::binsearch_insertpos(Base::hash_value_.direct_access(), cur_hash, size);

  if (pos >= size) 
    return false;
  
  const FlexibleStorage1D<Key>& cur_key_list = Base::key_[pos];
  const size_t sub_pos = Routines::find_unique(cur_key_list.direct_access(), key, cur_key_list.size());
  return (sub_pos < cur_key_list.size());
}

#endif