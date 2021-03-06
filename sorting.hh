/**** written by Thomas Schoenemann as a private person without employment, May 2013 ****/

#ifndef SORTING_HH
#define SORTING_HH

//NOTE: all routines are so far fixed to operator<

#ifndef MERGE_THRESH
#define MERGE_THRESH 256
#endif

#define QUICK_THRESH 8

static_assert(MERGE_THRESH >= 2);
static_assert(QUICK_THRESH >= 2);

#include "storage1D.hh"
#include "routines.hh"

//#include <cassert>
#include <algorithm> //for std::sort
#include <functional> //for std::less

//c++-11 provides this in <algorithm>
template<typename T, typename CmpOp = std::less<T> >
inline bool is_sorted(const T* data, const size_t nData) noexcept;

template<typename T, typename CmpOp = std::less<T> >
inline bool is_reverse_sorted(const T* data, const size_t nData) noexcept;

template<typename T>
inline bool is_unique_sorted(const T* data, const size_t nData) noexcept;

template<typename T>
inline bool is_unique_reverse_sorted(const T* data, const size_t nData) noexcept;

template<typename T, typename ST, typename CmpOp = std::less<T> >
bool is_index_sorted(const T* data, ST* indices, const ST nData) noexcept;

/***** plain sorting ****/

template<typename T, typename Swap = SwapOp<T>, typename CmpOp = std::less<T> >
inline void bubble_sort(T* data, const size_t nData) noexcept;

template<typename T, typename Swap = SwapOp<T>, typename CmpOp = std::less<T> >
inline void batch_bubble_sort(T* data, const size_t nData) noexcept;

template<typename T, typename Swap = SwapOp<T>, typename CmpOp = std::less<T> >
inline void bubble_sort_prepro(T* data, const size_t nData) noexcept;

template<typename T, typename CmpOp = std::less<T> > //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort(T* data, const size_t nData) noexcept;

template<typename T, typename CmpOp = std::less<T> > //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_withmoves(T* data, const size_t nData) noexcept;

//recursive, won't be inline!
template<typename T, typename Swap = SwapOp<T>, typename CmpOp = std::less<T> > //upshift can use c++ moves, no swap involved
void quick_sort(T* data, const size_t nData) noexcept;

//recursive, won't be inline!
template<typename T, typename Swap = SwapOp<T>, typename CmpOp = std::less<T>>
void merge_sort(T* data, const size_t nData) noexcept;

//interface to std::sort
template<typename T, typename CmpOp = std::less<T> >
inline void std_sort(T* data, const size_t nData) noexcept;

/***** key-value sort *****/

template<typename T1, typename T2, typename ValSwap = SwapOp<T2> >
inline void bubble_sort_key_value(T1* key, T2* value, const size_t nData) noexcept;

template<typename T1, typename T2, typename ValSwap = SwapOp<T2> >
inline void batch_bubble_sort_key_value(T1* key, T2* value, const size_t nData) noexcept;

template<typename T1, typename T2, typename ValSwap = SwapOp<T2> >
inline void bubble_sort_key_value_prepro(T1* key, T2* value, const size_t nData) noexcept;

template<typename T1, typename T2> //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_key_value(T1* key, T2* value, const size_t nData) noexcept;

template<typename T1, typename T2> //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_key_value_keymoves(T1* key, T2* value, const size_t nData) noexcept;

//recursive, won't be inline!
template<typename T1, typename T2, typename ValSwap = SwapOp<T2> >
void merge_sort_key_value(T1* key, T2* value, const size_t nData) noexcept;

/***** index sorting *****/

template<typename T, typename ST, typename CmpOp = std::less<T> >
inline void index_bubble_sort(const T* data, ST* indices, ST nData) noexcept;

template<typename T, typename ST, typename CmpOp = std::less<T> >
inline void index_batch_bubble_sort(const T* data, ST* indices, ST nData) noexcept;

//recursive, won't be inline!
template<typename T, typename ST, typename CmpOp = std::less<T> >
void index_merge_sort(const T* data, ST* indices, const ST nData) noexcept;

//recursive, won't be inline!
template<typename T, typename ST, typename CmpOp = std::less<T> >
void index_quick_sort(const T* data, ST* indices, const ST nData) noexcept;


/************ implementation  *************/

template<typename T, typename CmpOp>
inline bool is_sorted(const T* data, const size_t nData) noexcept
{
  static CmpOp less;
  
  //can we use AVX comparisons here for uint, int, float and double? 
  
  for (size_t i=1; i < nData; i++) {
    if (less(data[i],data[i-1]))
      return false;
  }

  return true;
}

template<typename T, typename CmpOp>
inline bool is_reverse_sorted(const T* data, const size_t nData) noexcept
{
  static CmpOp less;

  for (size_t i=1; i < nData; i++) {
    if (less(data[i-1],data[i]))
      return false;
  }

  return true;
}

template<typename T>
bool is_unique_sorted(const T* data, const size_t nData) noexcept
{  
  for (size_t i=1; i < nData; i++) {
    if (data[i] <= data[i-1])
      return false;
  }

  return true;
}

template<typename T>
inline bool is_unique_reverse_sorted(const T* data, const size_t nData) noexcept
{
  for (size_t i=1; i < nData; i++) {
    if (data[i-1] <= data[i])
      return false;
  }

  return true;
}

/**** plain ***/

template<typename T, typename Swap, typename CmpOp>
inline void bubble_sort(T* data, const size_t nData) noexcept
{
  const static Swap swapobj;
  const static CmpOp less;
  
  size_t last_flip = nData;

  while (last_flip > 1) {

    size_t flip = 0;

    for (size_t l=1; l < last_flip; l++) {

      if (less(data[l],data[l-1])) {
        swapobj(data[l],data[l-1]);
        flip = l;
      }
    }

    last_flip = flip;
  }
}

template<typename T, typename Swap, typename CmpOp>
inline void batch_bubble_sort(T* data, const size_t nData) noexcept
{ 
  const static CmpOp less;

  size_t last_flip = nData;

  while (last_flip > 1) {

    // std::cerr << "new iter, current data: ";
    // for (uint i=0; i < nData; i++)
      // std::cerr << data[i] << ", ";
    // std::cerr << std::endl;

    // std::cerr << "last flip: " << last_flip << std::endl;
    size_t flip = 0;

    for (size_t l=1; l < last_flip; l++) {

      if (less(data[l],data[l-1])) {

        //std::cerr << "found pair at " << l << std::endl;

        size_t end = l+1;
        for (; end < last_flip; end++) {

          if (!(less(data[end],data[end-1])))
            break;
        }

#ifdef USE_ROUTINES
        Routines::nontrivial_reverse<T,Swap>(data+l-1, end-l+1);
#else
        std::reverse(data+l-1, data+end);
#endif

        // std::cerr << "after reverse: ";
        // for (uint i=0; i < nData; i++)
          // std::cerr << data[i] << ", ";
        // std::cerr << std::endl;

        flip = end-1;
        l = end-2; //the last swapped index needs to be compared again (l is incremented by the loop)
      }
    }

    last_flip = flip;
  }
}

template<typename T, typename Swap, typename CmpOp>
inline void bubble_sort_prepro(T* data, const size_t nData) noexcept
{
  const static Swap swapobj;
  const static CmpOp less;

  const size_t thresh = 5;

  if (nData > thresh) {
    //try to ameliorate some bad cases
    for (uint i = nData-1; i >= thresh-1; i--) {
      if ( less(data[i],data[0]) )
        swapobj(data[i],data[0]);
    }
  }
}

template<typename T, typename CmpOp> //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort(T* data, const size_t nData) noexcept
{
  const static CmpOp less;
  
  for (size_t l=1; l < nData; l++) {
    const T curdat = data[l];
    size_t inspos = l;
    while (inspos > 0 && less(curdat,data[inspos-1]))
      inspos--;

    if (inspos != l) {
      Routines::upshift_array(data, inspos, l, 1);
      //for (size_t ll = l; ll > inspos; ll--)
      //  data[ll] = data[ll-1];
      data[inspos] = curdat;
    }
  }
}

template<typename T, typename CmpOp> //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_withmoves(T* data, const size_t nData) noexcept
{
  const static CmpOp less;
  
  for (size_t l=1; l < nData; l++) {
    const T& curdat = data[l];
    size_t inspos = l;
    while (inspos > 0 && less(curdat,data[inspos-1]))
      inspos--;

    if (inspos != l) {
      T temp = std::move(data[l]);
      Routines::upshift_array(data, inspos, l, 1);
      //for (size_t ll = l; ll > inspos; ll--)
      //  data[ll] = data[ll-1];
      data[inspos] = temp;
    }
  }
}

//recursive, won't be inline!
template<typename T, typename Swap, typename CmpOp> //upshift can use c++ moves, no swap involved
void quick_sort(T* data, const size_t nData) noexcept
{
  const static Swap swapobj;
  const static CmpOp less;

  if (nData <= 1)
    return;
  if (nData == 2) {
    if (data[1] < data[0])
      swapobj(data[0],data[1]);
  }
  else if (nData <= QUICK_THRESH)
    bubble_sort<T,Swap,CmpOp>(data,nData);
  else {

  //std::cerr << "quick_sort(" << nData << ")" << std::endl;

    using KeyRefType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;

    KeyRefType pivot = data[nData-1];

    size_t i=0;
    size_t j=nData-2; //not the perfect type (in view of ST)

    while (true) {

      while(i < int(nData-1) && !(less(pivot,data[i])) ) // data[i] <= pivot)
        i++;

      while(j > i && !(less(data[j],pivot)) ) // data[j] >= pivot)
        j--;

      if (i >= j)
        break;

      std::swap(data[i],data[j]);
      i++;
      j--;
    }

    if (i != nData-1)
      swapobj(data[i],data[nData-1]);
  
    //recursive calls
    if (i > 1) { //no need to sort a single entry
      quick_sort<T,Swap,CmpOp>(data,i);
    }
    if (i < nData-2) {
      quick_sort<T,Swap,CmpOp>(data+i+1,nData-i-1);
    }
  }
  
  //bool issorted = is_sorted<T,CmpOp>(data, nData);
  //assert(issorted);
}

template<typename T, typename SwapOp, typename CmpOp>
void merge_sort(T* data, const size_t nData) noexcept
{  
  //merge sort is meant for large arrays => check this for speed
  if (is_sorted<T,CmpOp>(data, nData))
    return;

  if (nData <= MERGE_THRESH) 
    bubble_sort<T,SwapOp,CmpOp>(data, nData);
  else {
 
    const static CmpOp less;
 
    using RefType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;

    const size_t half = nData / 2;
    const size_t nData2 = nData-half;

    T* aux_data = new T[nData];

    Makros::unified_move_assign(aux_data, data, nData);
    
    //Makros::unified_move_assign(aux_data, data, half);
    //Makros::unified_move_assign(aux_data+half, data+half, nData2);

    merge_sort<T,SwapOp,CmpOp>(aux_data,half);
    merge_sort<T,SwapOp,CmpOp>(aux_data+half,nData2);
    
    //now merge the lists
    const T* data1 = aux_data;
    const T* data2 = aux_data+half;

    size_t k1=0;
    size_t k2=0;
    size_t k=0;

    while (true) {

      RefType d1 = data1[k1];
      RefType d2 = data2[k2];
            
      if (!less(d2,d1)) {
        data[k] = std::move(d1);
        k1++;
        if (k1 >= half)
          break;
      }
      else {
        data[k] = std::move(d2);
        k2++;
        if (k2 >= nData2)
          break;
      }

      k++;
    }

    //copy the rest
    k++;
    Makros::unified_move_assign(data+k, data1+k1, half-k1);
    Makros::unified_move_assign(data+k, data2+k2, nData2-k2);

    delete[] aux_data;
  }
}

//interface to std::sort
template<typename T, typename CmpOp>
inline void std_sort(T* data, const size_t nData) noexcept
{
  std::sort<T,CmpOp>(data, data+nData);
}

/**** key-value ****/

template<typename T1, typename T2, typename ValSwap>
inline void bubble_sort_key_value(T1* key, T2* value, const size_t nData) noexcept
{
  const static ValSwap valswapobj;
  
  size_t last_flip = nData;

  while (last_flip > 1) {

    size_t flip = 0;

    for (size_t j=1; j < last_flip; j++) {

      if (key[j] < key[j-1]) {

        std::swap(key[j],key[j-1]);
        valswapobj(value[j],value[j-1]);
        flip = j;
      }
    }

    last_flip = flip;
  }
}

template<typename T1, typename T2, typename ValSwap>
inline void batch_bubble_sort_key_value(T1* key, T2* value, const size_t nData) noexcept
{
  size_t last_flip = nData;

  while (last_flip > 1) {

    //std::cerr << "*** keys: ";
    //for (uint k = 0; k < nData; k++)
    //  std::cerr << key[k] << " ";
    //std::cerr << std::endl;
    //std::cerr << "last flip: " << last_flip << std::endl;

    size_t flip = 0;

    for (size_t l=1; l < last_flip; l++) {

      //std::cerr << "l: " << l << std::endl;

      if (key[l] < key[l-1]) {

        size_t end = l+1;
        for (; end < last_flip; end++) {
          if (!(key[end] < key[end-1]))
            break;
        }

        //std::cerr << "reversing from " << (l-1) << " to excl. " << end << std::endl;
#ifdef USE_ROUTINES
        Routines::nontrivial_reverse(key+l-1, end-l+1);
        Routines::nontrivial_reverse<T2,ValSwap>(value+l-1, end-l+1);
#else
        std::reverse(key+l-1, key+end);
        std::reverse(value+l-1, value+end);
#endif        
        flip = end-1;
        l = end-2; //the last swapped index needs to be compared again (l is incremented by the loop)
      }
    }

    //std::cerr << "*** keys after: ";
    //for (uint k = 0; k < nData; k++)
    //  std::cerr << key[k] << " ";
    //std::cerr << std::endl;

    last_flip = flip;
  }
}

template<typename T1, typename T2, typename ValSwap>
inline void bubble_sort_key_value_prepro(T1* key, T2* value, const size_t nData) noexcept
{
  const static ValSwap valswapobj;

  const size_t thresh = 5;

  if (nData > thresh) {
    //try to ameliorate some bad cases
    for (uint i = nData-1; i >= thresh-1; i--) {
      if (key[i] < key[0]) {
        std::swap(key[i],key[0]);
        valswapobj(value[i],value[0]);
      }
    }
  }
}

template<typename T1, typename T2>
inline void shift_bubble_sort_key_value(T1* key, T2* value, const size_t nData) noexcept
{  
  for (size_t j=1; j < nData; j++) {
    const T1 curkey = key[j];
    size_t inspos = j;
    while (inspos > 0 && curkey < key[inspos-1])
      inspos--;

    if (inspos != j) {
      const T2 curval = std::move(value[j]);
      Routines::upshift_array(key, inspos, j, 1);
      Routines::upshift_array(value, inspos, j, 1);
      // for (size_t jj = j; jj > inspos; jj--) {
        // key[jj] = key[jj-1];
        // value[jj] = value[jj-1];
      // }
      key[inspos] = curkey;
      value[inspos] = curval;
    }
  }
}

template<typename T1, typename T2> //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_key_value_keymoves(T1* key, T2* value, const size_t nData) noexcept
{
  for (size_t j=1; j < nData; j++) {
    const T1 &curkey = key[j];
    size_t inspos = j;
    while (inspos > 0 && curkey < key[inspos-1])
      inspos--;

    if (inspos != j) {
      const T1 temp = std::move(key[j]);
      const T2 curval = std::move(value[j]);
      Routines::upshift_array(key, inspos, j, 1);
      Routines::upshift_array(value, inspos, j, 1);
      // for (size_t jj = j; jj > inspos; jj--) {
        // key[jj] = key[jj-1];
        // value[jj] = value[jj-1];
      // }
      key[inspos] = temp;
      value[inspos] = curval;
    }
  }  
}

//recursive, won't be inline!
template<typename T1, typename T2, typename ValSwap>
void merge_sort_key_value(T1* key, T2* value, const size_t nData) noexcept
{
  //merge sort is meant for large arrays => check this for speed
  if (is_sorted(key, nData))
    return;
  
  if (nData <= MERGE_THRESH)
    bubble_sort_key_value<T1,T2,ValSwap>(key, value, nData);
  else {

    using KeyRefType = typename std::conditional<std::is_fundamental<T1>::value || std::is_pointer<T1>::value, const T1, const T1&>::type;

    const size_t half = nData / 2;
    const size_t nData2 = nData-half;

    T1* aux_key = new T1[nData];
    T2* aux_value = new T2[nData];

    Makros::unified_move_assign(aux_key, key, nData);
    Makros::unified_move_assign(aux_value, value, nData);

    //Makros::unified_move_assign(aux_key, key, half);
    //Makros::unified_move_assign(aux_value, value, half);

    merge_sort_key_value(aux_key,aux_value,half);

    //Makros::unified_move_assign(aux_key+half, key+half, nData2);
    //Makros::unified_move_assign(aux_value+half, value+half, nData2);

    merge_sort_key_value(aux_key+half,aux_value+half,nData2);

    //now merge the lists

    const T1* key1 = aux_key;
    const T1* key2 = aux_key+half;
    const T2* value1 = aux_value;
    const T2* value2 = aux_value+half;

    size_t k1=0;
    size_t k2=0;
    size_t k=0;

    while (true) {

      const KeyRefType d1 = key1[k1];
      const KeyRefType d2 = key2[k2];
      if (d1 <= d2) {
        key[k] = std::move(d1);
        value[k] = std::move(value1[k1]);
        k1++;
        if (k1 >= half)
          break;
      }
      else {
        key[k] = std::move(d2);
        value[k] = std::move(value2[k2]);
        k2++;
        if (k2 >= nData2)
          break;
      }

      k++;
    }

    //copy the rest
    k++;
    Makros::unified_move_assign(key+k, key1+k1, half-k1);
    Makros::unified_move_assign(value+k, value1+k1, half-k1);
    Makros::unified_move_assign(key+k, key2+k2, nData2-k2);
    Makros::unified_move_assign(value+k, value2+k2, nData2-k2);

    delete[] aux_key;
    delete[] aux_value;
  }
}

/**** index ****/

template<typename T, typename ST, typename CmpOp>
bool is_index_sorted(const T* data, ST* indices, const ST nData) noexcept
{
  const static CmpOp less;
  
  for (ST i = 1; i < nData; i++) {
    
    if (less(data[indices[i]],data[indices[i-1]]))
      return false;
  }
  
  return true;
}

template<typename T, typename ST, typename CmpOp = std::less<T> >
inline void aux_index_bubble_sort(const T* data, ST* indices, const ST nData) noexcept;

template<typename T, typename ST, typename CmpOp>
inline void index_bubble_sort(const T* data, ST* indices, const ST nData) noexcept
{
  for (uint i=0; i < nData; i++)
    indices[i] = i;

  aux_index_bubble_sort<T,ST,CmpOp>(data, indices, nData);
  
#ifndef NDEBUG    
  bool is_sorted = is_index_sorted<T,ST,CmpOp>(data, indices, nData);
  assert(is_sorted);
#endif  
}

template<typename T, typename ST, typename CmpOp>
inline void aux_index_bubble_sort(const T* data, ST* indices, const ST nData) noexcept
{
  const static CmpOp less;

  ST last_flip = nData;

  while (last_flip > 1) {

    ST flip = 0;

    for (uint l=1; l < last_flip; l++) {

      if (less(data[indices[l]],data[indices[l-1]])) {
        std::swap(indices[l],indices[l-1]);
        flip = l;
      }
    }
    
    last_flip = flip;
  }
}

template<typename T, typename ST, typename CmpOp>
inline void index_batch_bubble_sort(const T* data, ST* indices, ST nData) noexcept
{
  const static CmpOp less;

  for (uint i=0; i < nData; i++)
    indices[i] = i;

  size_t last_flip = nData;

  while (last_flip > 1) {

    size_t flip = 0;

    for (size_t l=1; l < last_flip; l++) {

      if (less(data[indices[l]],data[indices[l-1]])) {

        size_t end = l+1;
        for (; end < last_flip; end++) {

          if (!(less(data[indices[end]],data[indices[end-1]])))
            break;
        }

        std::reverse(indices+l-1, indices+end);

        flip = l;
        l = end-2; //the last swapped index needs to be compared again (l is incremented by the loop)
      }
    }

    last_flip = flip;
  }

#ifndef NDEBUG    
  bool is_sorted = is_index_sorted<T,ST,CmpOp>(data, indices, nData);
  assert(is_sorted);
#endif  
}

//recursive, won't be inline!
template<typename T, typename ST, typename CmpOp = std::less<T> >
void aux_index_quick_sort(const T* data, ST* indices, const ST nData) noexcept;

//recursive, won't be inline!
template<typename T, typename ST, typename CmpOp>
void aux_index_merge_sort(const T* data, ST* indices, const ST nData) noexcept
{
  //std::cerr << "nData: "<< nData << std::endl;

  // if (nData <= 12 /*8*/) {

  //   //bubble sort
  //   for (uint k=0; k < nData-1; k++) {

  //     //std::cerr << "k: " << k << std::endl;

  //     for (uint l=0; l < nData-1-k; l++) {
  // 	//std::cerr << "l: " << l << std::endl;

  // 	//std::cerr << "index1: " << indices[l] << std::endl;
  // 	//std::cerr << "index2: " << indices[l+1] << std::endl;

  // 	if (data[indices[l]] > data[indices[l+1]])
  // 	  std::swap(indices[l],indices[l+1]);
  //     }
  //   }
  // }

  //merge sort is meant for large arrays => check this for speed
  if (is_sorted(data, nData))
    return;

  if (nData <= MERGE_THRESH) {
    aux_index_quick_sort<T,ST,CmpOp>(data,indices,nData);
  }
  else {

    const static CmpOp less;

    const ST half = nData / 2;
    const ST nData2 = nData-half;

    ST* aux_indices = new ST[nData];

    Makros::unified_move_assign(aux_indices, indices, nData);

    //Makros::unified_move_assign(aux_indices, indices, half);

    aux_index_merge_sort<T,ST,CmpOp>(data,aux_indices,half);

    //Makros::unified_move_assign(aux_indices+half, indices+half, nData2);

    aux_index_merge_sort<T,ST,CmpOp>(data,aux_indices+half,nData2);

    const ST* index1 = aux_indices;
    const ST* index2 = aux_indices+half;

    ST k1=0;
    ST k2=0;
    ST k=0;

    //while(k1 < half && k2 < nData2) {
    while (true) {

      const ST idx1 = index1[k1];
      const ST idx2 = index2[k2];
      if (!(less(data[idx2],data[idx1]))) {
        indices[k] = idx1;
        k1++;
        if (k1 >= half)
          break;
      }
      else {
        indices[k] = idx2;
        k2++;
        if (k2 >= nData2)
          break;
      }

      k++;
    }

    //copy the rest
    k++;
    Makros::unified_move_assign(indices+k, index1+k1, half-k1);
    Makros::unified_move_assign(indices+k,index2+k2,nData2-k2);

    delete[] aux_indices;
  }
}

//recursive, won't be inline!
//this is still fixed to memcpy
template<typename T, typename ST, typename CmpOp>
void aux_index_merge_sort_4split(const T* data, ST* indices, const ST nData) noexcept
{
  //std::cerr << "nData: "<< nData << std::endl;

  if (nData <= 32) {
    aux_index_quick_sort<T,ST,CmpOp>(data,indices,nData);
  }
  else {

    const static CmpOp less;

    uint nData1to3 = nData/4;
    uint nData4 = nData-3*nData1to3;

    const ST k_limit[4] = {nData1to3,nData1to3,nData1to3,nData4};

    //TODO: use unified_move_assign and a single array
    //ST* aux_indices = new ST[nData];
    //Makros::unified_move_assign(aux_indices, indices, nData);
    //Storage1D<ST*> data[4] = {aux_indices,aux_indices+nData1to3,aux_indices+2*nData1to3,aux_indices3*nData1to3};

    Storage1D<ST> index[4];

    index[0].resize(nData1to3);
    memcpy(index[0].direct_access(),indices,nData1to3*sizeof(ST));

    aux_index_merge_sort<T,ST,CmpOp>(data,index[0].direct_access(),nData1to3);

    index[1].resize(nData1to3);
    memcpy(index[1].direct_access(),indices+nData1to3,nData1to3*sizeof(ST));

    aux_index_merge_sort<T,ST,CmpOp>(data,index[1].direct_access(),nData1to3);

    index[2].resize(nData1to3);
    memcpy(index[2].direct_access(),indices+2*nData1to3,nData1to3*sizeof(ST));

    aux_index_merge_sort<T,ST,CmpOp>(data,index[2].direct_access(),nData1to3);

    index[3].resize(nData4);
    memcpy(index[3].direct_access(),indices+3*nData1to3,nData4*sizeof(ST));

    aux_index_merge_sort<T,ST,CmpOp>(data,index[3].direct_access(),nData4);

    ST k_idx[4] = {0,0,0,0};

    uint k=0;
    while (k < nData) {

      uint arg_best = MAX_UINT;
      uint idx;
      T best;
      for (uint kk=0; kk < 4; kk++) {

        const uint cur_idx = k_idx[kk];

        if (cur_idx < k_limit[kk]) {

          const uint real_idx = index[kk][cur_idx];

          const T cur_val = data[real_idx];

          if (arg_best == MAX_UINT || less(cur_val,best)) {
            best = cur_val;
            arg_best=kk;
            idx = real_idx;
          }
        }
      }

      indices[k] = idx;
      k_idx[arg_best]++;
      k++;
    }
  }
}

//recursive, won't be inline!
template<typename T, typename ST, typename CmpOp>
void index_merge_sort(const T* data, ST* indices, const ST nData) noexcept
{
  //std::cerr << "sorting " << nData << " entries" << std::endl;

  for (uint i=0; i < nData; i++)
    indices[i] = i;

  aux_index_merge_sort<T,ST,CmpOp>(data,indices,nData);
  
#ifndef NDEBUG    
  bool is_sorted = is_index_sorted<T,ST,CmpOp>(data, indices, nData);
  assert(is_sorted);
#endif  
}


//recursive, won't be inline!
template<typename T, typename ST, typename CmpOp>
void aux_index_quick_sort(const T* data, ST* indices, const ST nData) noexcept
{
  //std::cerr << "nData: " << nData << std::endl;

  if (nData <= QUICK_THRESH) 
  {
    aux_index_bubble_sort<T,ST,CmpOp>(data, indices, nData);
  }
  else {

    using KeyRefType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;

    const static CmpOp less;

    KeyRefType pivot = data[indices[nData-1]];

    ST i=0;
    ST j=nData-2; 

    // if (nData <= 75) {
    //   std::cerr << "sequence: ";
    //   for (uint k=0; k < nData; k++)
    // 	std::cerr << data[indices[k]] << " ";
    //   std::cerr << std::endl;
    // }

    while (true) {

      while(i < (nData-1) && !(less(pivot,data[indices[i]])) ) // data[indices[i]] <= pivot)
        i++;

      while(j > i && !(less(data[indices[j]],pivot)) ) // data[indices[j]] >= pivot)
        j--;

      if (i >= j)
        break;

      std::swap(indices[i],indices[j]);
      i++;
      j--;
    }

    if (i != nData-1)
      std::swap(indices[i],indices[nData-1]);


    // if (nData <= 75) {
    //   std::cerr << "sequence after shifting: ";
    //   for (uint k=0; k < nData; k++)
    // 	std::cerr << data[indices[k]] << " ";
    //   std::cerr << std::endl;
    // }


    //std::cerr << "i: " << i << ", j: " << j << std::endl;

    //recursive calls

    if (i > 1) { //no need to sort a single entry
      aux_index_quick_sort<T,ST,CmpOp>(data,indices,i);
    }
    if (i < nData-2) {
      aux_index_quick_sort<T,ST,CmpOp>(data,indices+i+1,nData-i-1);
    }
  }

}

//recursive, won't be inline!
template<typename T, typename ST, typename CmpOp>
void index_quick_sort(const T* data, ST* indices, const ST nData) noexcept
{
  //std::cerr << "sorting " << nData << " entries" << std::endl;

  for (uint i=0; i < nData; i++)
    indices[i] = i;

  aux_index_quick_sort<T,ST,CmpOp>(data,indices,nData);

#ifndef NDEBUG    
  bool is_sorted = is_index_sorted<T,ST,CmpOp>(data, indices, nData);
  assert(is_sorted);
#endif  
}


#endif
