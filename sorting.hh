/**** written by Thomas Schoenemann as a private person without employment, May 2013 ****/

#ifndef SORTING_HH
#define SORTING_HH

#ifndef MERGE_THRESH
#define MERGE_THRESH 256
#endif

#ifndef QUICK_THRESH
#define QUICK_THRESH 16
#endif

#define USE_ROUTINES


#include "storage1D.hh"
#include "routines.hh"
#include "flexible_heap_priority_queue.hh"

#include <algorithm>

//c++-11 provides this in <algorithm>
template<typename T, typename Less = std::less<T> >
bool is_sorted(const T* data, const size_t nData) noexcept;

template<typename T, typename Less = std::less<T> >
bool is_reverse_sorted(const T* data, const size_t nData) noexcept;

template<typename T, typename Less = std::less<T>, typename Equal = std::equal_to<T>  >
bool is_unique_sorted(const T* data, const size_t nData) noexcept;

template<typename T, typename Less = std::less<T>, typename Equal = std::equal_to<T> >
bool is_unique_reverse_sorted(const T* data, const size_t nData) noexcept;

template<typename T, typename ST, typename Less = std::less<T> >
bool is_index_sorted(T* data, ST* indices, const size_t nData) noexcept;


/***** plain sorting ****/

template<typename T, typename Swap = SwapOp<T>, typename Less = std::less<T> >
inline void bubble_sort(T* data, const size_t nData) noexcept;

template<typename T, typename Swap = SwapOp<T>, typename Less = std::less<T> >
inline void batch_bubble_sort(T* data, const size_t nData) noexcept;

template<typename T, typename Swap = SwapOp<T>, typename Less = std::less<T> >
inline void bubble_sort_prepro(T* data, const size_t nData) noexcept;

template<typename T, typename Swap = SwapOp<T>, typename Less = std::less<T> >
inline void shaker_sort(T* data, const size_t nData) noexcept;

template<typename T, typename Less = std::less<T> > //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort(T* data, const size_t nData) noexcept;

template<typename T, typename Less = std::less<T> > //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_withmoves(T* data, const size_t nData) noexcept;

template<typename T, typename Less = std::less<T>, typename Equal = std::equal_to<T>, bool stable = true >
inline void insertion_sort(T* data, const size_t nData) noexcept;

//fast, but not stable and not in-place (actually tree sort)
template<typename T, typename Less = std::less<T> >
inline void selection_sort(T* data, const size_t nData) noexcept;

//fast, but not stable and not in-place (actually tree sort)
template<typename T, typename Less = std::less<T> >
inline void selection_sort_moves(T* data, const size_t nData) noexcept;


//fast and in-place, but not stable
template<typename T, typename Less = std::less<T> >
inline void heap_sort(T* data, const size_t nData) noexcept;

//fast and in-place, but not stable
template<typename T, typename Less = std::less<T> >
inline void heap_sort_moves(T* data, const size_t nData) noexcept;

//recursive, won't be inline
//fallback_algo: 0 = insertion_sort, 1 = bubble sort, 2 = batch bubble sort
template<typename T, typename Swap = SwapOp<T>, typename Less = std::less<T>, typename Equal = std::equal_to<T>, int fallback_algo = 0 >
void quick_sort(T* data, const size_t nData) noexcept;

//recursive, won't be inline
template<typename T, typename Swap = SwapOp<T>, typename Less = std::less<T>, typename Equal = std::equal_to<T> >
void quick_sort_split_equals(T* data, const size_t nData) noexcept;

//recursive, won't be inline
template<typename T, typename Swap = SwapOp<T>, typename Less = std::less<T>, typename Equal = std::equal_to<T> >
void merge_sort(T* data, const size_t nData) noexcept;

//recursive, won't be inline
template<typename T, typename Swap = SwapOp<T>, typename Less = std::less<T>, typename Equal = std::equal_to<T> >
void qmc_sort(T* data, const size_t nData) noexcept;

//interface to std::sort
template<typename T>
void std_sort(T* data, const size_t nData) noexcept;

/***** key-value sort *****/

template<typename TK, typename TV, typename ValSwap = SwapOp<TV>, typename Less = std::less<TK> >
inline void bubble_sort_key_value(TK* key, TV* value, const size_t nData) noexcept;

template<typename TK, typename TV, typename ValSwap = SwapOp<TV>, typename Less = std::less<TK> >
inline void batch_bubble_sort_key_value(TK* key, TV* value, const size_t nData) noexcept;

template<typename TK, typename TV, typename ValSwap = SwapOp<TV>, typename Less = std::less<TK> >
inline void bubble_sort_key_value_prepro(TK* key, TV* value, const size_t nData) noexcept;

template<typename TK, typename TV, typename Less = std::less<TK> > //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_key_value(TK* key, TV* value, const size_t nData) noexcept;

template<typename TK, typename TV, typename Less = std::less<TK> > //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_key_value_keymoves(TK* key, TV* value, const size_t nData) noexcept;

template<typename TK, typename TV, typename Less = std::less<TK>, typename Equal = std::equal_to<TK>, bool stable = true>
inline void insertion_sort_key_value(TK* key, TV* value, const size_t nData) noexcept;

//not stable and not in-place
template<typename TK, typename TV, typename Less = std::less<TK> >
inline void selection_sort_key_value(TK* key, TV* value, const size_t nData) noexcept;

//in-place, but not stable
template<typename TK, typename TV, typename Less = std::less<TK> >
inline void heap_sort_key_value(TK* key, TV* value, const size_t nData) noexcept;


//recursive, won't be inline
template<typename TK, typename TV, typename Swap = SwapOp<TV>, typename Less = std::less<TK>, typename Equal = std::equal_to<TK> >
void quick_sort_key_value(TK* data, TV* value, const size_t nData) noexcept;

//recursive, won't be inline
template<typename TK, typename TV, typename Swap = SwapOp<TV>, typename Less = std::less<TK>, typename Equal = std::equal_to<TK> >
void merge_sort_key_value(TK* key, TV* value, const size_t nData) noexcept;

//recursive, won't be inline
template<typename TK, typename TV, typename ValSwap = SwapOp<TV>, typename Less = std::less<TK>, typename Equal = std::equal_to<TK> >
void qmc_sort_key_value(TK* data, TV* value, const size_t nData) noexcept;

/***** index sorting *****/

template<typename T, typename ST>
inline void index_bubble_sort(const T* data, ST* indices, ST nData) noexcept;

template<typename T, typename ST>
inline void index_batch_bubble_sort(const T* data, ST* indices, ST nData) noexcept;

template<typename T, typename ST, bool stable = true>
inline void index_insertion_sort(T* data, ST* indices, const ST nData) noexcept;

template<typename T, typename ST>
inline void index_selection_sort(const T* data, ST* indices, ST nData) noexcept;

//recursive, won't be inline
template<typename T, typename ST>
void index_quick_sort(const T* data, ST* indices, const ST nData) noexcept;

//recursive, won't be inline
template<typename T, typename ST>
void index_merge_sort(const T* data, ST* indices, const ST nData) noexcept;



/************ implementation  *************/

template<typename T, typename Less>
bool is_sorted(const T* data, const size_t nData) noexcept
{
  //can we use AVX comparisons here for uint, int, float and double?

  const static Less less;

  for (size_t i=1; i < nData; i++) {
    if (less(data[i],data[i-1]))
      return false;
  }

  return true;
}

template<typename T, typename Less>
bool is_reverse_sorted(const T* data, const size_t nData) noexcept
{
  const static Less less;

  for (size_t i=1; i < nData; i++) {
    if (less(data[i-1],data[i]))
      return false;
  }

  return true;
}

template<typename T, typename Less, typename Equal>
bool is_unique_sorted(const T* data, const size_t nData) noexcept
{
  const static Less less;
  const static Equal equal;

  for (size_t i=1; i < nData; i++) {
    if (less(data[i],data[i-1]) || equal(data[i],data[i-1]))
      return false;
  }

  return true;
}

template<typename T, typename Less, typename Equal>
bool is_unique_reverse_sorted(const T* data, const size_t nData) noexcept
{
  const static Less less;
  const static Equal equal;

  for (size_t i=1; i < nData; i++) {
    if (less(data[i-1],data[i]) || equal(data[i-1],data[i]))
      return false;
  }

  return true;
}

template<typename T, typename ST, typename Less = std::less<T> >
bool is_index_sorted(T* data, ST* indices, const size_t nData) noexcept
{
  const static Less less;
  for (ST i = 1; i < nData; i++) {
    if (less(data[indices[i]],data[indices[i-1]]))
      return false;
  }
  return true;
}



/**** plain ***/

template<typename T, typename Swap, typename Less>
inline void bubble_sort(T* data, const size_t nData) noexcept
{
  const static Swap swapobj;
  const static Less less;

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

template<typename T, typename Swap, typename Less>
inline void batch_bubble_sort(T* data, const size_t nData) noexcept
{
  const static Less less;

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

          if (!less(data[end],data[end-1]))
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

template<typename T, typename Swap>
inline void bubble_sort_prepro(T* data, const size_t nData) noexcept
{
  const static Swap swapobj;

  const size_t thresh = 5;

  if (nData > thresh) {
    //try to ameliorate some bad cases
    for (uint i = nData-1; i >= thresh-1; i--) {
      if ( data[i] < data[0] )
        swapobj(data[i],data[0]);
    }
  }
}

template<typename T, typename Less> //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort(T* data, const size_t nData) noexcept
{
  const static Less less;

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

template<typename T, typename Swap = SwapOp<T>, typename Less = std::less<T> >
inline void shaker_sort(T* data, const size_t nData) noexcept
{
  const static Swap swapobj;
  const static Less less;

  size_t last_flip = nData;
  size_t first_flip = 1;

  while (last_flip > 1) {

    size_t flip = 0;

    for (size_t l=first_flip; l < last_flip; l++) {

      if (less(data[l],data[l-1])) {
        swapobj(data[l],data[l-1]);
        flip = l;
      }
    }

    last_flip = flip;
    flip = 1;
    for (size_t l = last_flip; l >= first_flip; l--) {
      if (less(data[l],data[l-1])) {
        swapobj(data[l],data[l-1]);
        flip = l;
      }
    }

    first_flip = flip;
  }
}

template<typename T, typename Less> //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_withmoves(T* data, const size_t nData) noexcept
{
  const static Less less;
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

template<typename T, typename Less, typename Equal, bool stable>
inline void insertion_sort(T* data, const size_t nData) noexcept
{
  if (nData <= 1)
    return;

  const static Less less;
  const static Equal equal;

  if (less(data[1],data[0]))
    std::swap(data[0],data[1]);

  for (size_t i = 2; i < nData; i++) {

    size_t k = Routines::binsearch_insertpos<T,Less,Equal>(data, data[i], i);
    if (stable) {
      while (k < i && equal(data[k],data[i]))
        k++;
    }

    if (k < i) {
      T temp = std::move(data[i]);
      Routines::upshift_array(data, k, i, 1);
      data[k] = std::move(temp);
    }
  }
}

template<typename T, typename Less>
inline void selection_sort(T* data, const size_t nData) noexcept
{
  if (nData < 2)
	return; //nothing to do
	
  FlexibleHeapPriorityQueue<T,Less> heap_queue(nData);
  
  //std::cerr << "init heap" << std::endl;
  heap_queue.insert_first_pair(std::move(data[0]),std::move(data[1]));
  for (size_t idx = 2; idx < nData; idx++)
    heap_queue.insert(data[idx]);
 
  //std::cerr << "proceed" << std::endl;
  for (size_t idx = 0; idx < nData-2; idx++)
    data[idx] = heap_queue.extract_min();
  heap_queue.extract_last_pair(data[nData-2],data[nData-1]);
}

template<typename T, typename Less>
inline void selection_sort_moves(T* data, const size_t nData) noexcept
{
  if (nData < 2)
	return; //nothing to do
	
  FlexibleHeapPriorityQueue<T,Less> heap_queue(nData);
  
  //std::cerr << "init heap" << std::endl;
  heap_queue.insert_first_pair(std::move(data[0]),std::move(data[1]));
  for (size_t idx = 2; idx < nData; idx++)
    heap_queue.move_insert(std::move(data[idx]));	  
 
  //std::cerr << "proceed" << std::endl;
  for (size_t idx = 0; idx < nData-2; idx++)
    data[idx] = heap_queue.move_extract_min();
  heap_queue.extract_last_pair(data[nData-2],data[nData-1]);
}

//fast and in-place, but not stable
template<typename T, typename Less>
inline void heap_sort(T* data, const size_t nData) noexcept
{
  if (nData < 2)
	return; //nothing to be done

  //NOTE: heap is sorted upwards!

  T* heap = data-1;
  const static Less less;
  
  if (less(heap[1],heap[2]))
    std::swap(heap[1],heap[2]);

  //forward - insert nodes into the heap
  for(size_t idx = 3; idx <= nData; idx++) {
    //insert heap[idx] into the heap. it is already at the last position
    const T key = std::move(heap[idx]);
	
	size_t i = idx;
	
    while (i > 1) {

      const size_t parent = i >> 1; //i / 2;
      if (less(heap[parent],key)) {
        heap[i] = std::move(heap[parent]);
	    i = parent;
      }
      else
        break;
    } 	
  
    heap[i] = key;
  }	  

  //backward - dequeue nodes from heap
  for (size_t idx = nData; idx > 2; idx--)
  {	
	const T key = std::move(heap[idx]);
	heap[idx] = std::move(heap[1]);
	
	size_t i = 1;
	
    while (true) {

      const size_t j1 = i << 1; //2*i;

      if (j1 >= idx)
        break;

      size_t better_j = j1;

	  if (j1 == idx-1) {
	    //only one child exists
	  }
      else {
        //both children exist

        const size_t j2 = j1 | 1; // j1+1;

        if (less(heap[j1],heap[j2]))
          better_j = j2;
	  }

      if (less(key,heap[better_j])) {
        heap[i] = std::move(heap[better_j]); 
	    i = better_j;
      }
      else
        break;
	}
	
    heap[i] = key;
  }

  std::swap(heap[1],heap[2]);  
}

//fast and in-place, but not stable
template<typename T, typename Less>
inline void heap_sort_moves(T* data, const size_t nData) noexcept
{
  if (nData < 2)
	return; //nothing to be done

  //NOTE: heap is sorted upwards!

  T* heap = data-1;
  const static Less less;
  
  if (less(heap[1],heap[2]))
    std::swap(heap[1],heap[2]);

  //forward - insert nodes into the heap
  for(size_t idx = 3; idx <= nData; idx++) {
    //insert heap[idx] into the heap. it is already at the last position
	T key = std::move(heap[idx]);
	const T& ckey = key;
	
	size_t i = idx;
	
    while (i > 1) {

      const size_t parent = i >> 1; //i / 2;
      if (less(heap[parent],ckey)) {
        heap[i] = std::move(heap[parent]);
	    i = parent;
      }
      else
        break;
    } 	
  
    heap[i] = std::move(key);  	  
  }
  
  //backward - dequeue nodes from heap
  for (size_t idx = nData; idx > 2; idx--)
  {	
	T key = std::move(heap[idx]);
	const T& ckey = key;
	heap[idx] = std::move(heap[1]);
	
	size_t i = 1;
	
    while (true) {

      const size_t j1 = i << 1; //2*i;

      if (j1 >= idx)
        break;

      size_t better_j = j1;

	  if (j1 == idx-1) {
	    //only one child exists
	  }
      else {
        //both children exist

        const size_t j2 = j1 | 1; //j1+1;

        if (less(heap[j1],heap[j2]))
          better_j = j2;
	  }

      if (less(ckey,heap[better_j])) {
        heap[i] = std::move(heap[better_j]); 
	    i = better_j;
      }
      else
        break;
	}
	
    heap[i] = std::move(key);
  }

  std::swap(heap[1],heap[2]);  
}



//fallback_algo: 0 = insertion_sort, 1 = bubble sort, 2 = batch bubble sort
template<typename T, typename Swap, typename Less, typename Equal, int fallback_algo>
void quick_sort(T* data, const size_t nData) noexcept
{
  if (nData <= QUICK_THRESH) {
    if (fallback_algo == 0)
      insertion_sort<T,Less,Equal,true>(data,nData);
    else if (fallback_algo == 1)
      bubble_sort<T,Swap,Less>(data,nData);
    else if (fallback_algo == 2)
      batch_bubble_sort<T,Swap,Less>(data,nData);
  }
  else {

    const static Less less;
    const static Swap swapobj;
    using PassType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;
    PassType pivot = data[nData-1];

    size_t i=0;
    size_t j=nData-2; //not the perfect type (in view of ST)

    while(true) {

      //while((i < j) && (pivot < data[j]))
      while((i < j) && less(pivot,data[j]))
        j--;

      //while((i < j) && (data[i] <= pivot))
      while((i < j) && !less(pivot,data[i]))
        i++;

      if (i >= j)
        break;

      assert(less(data[j],data[i]));
      swapobj(data[i],data[j]);
#if 0
      i++;
      j--;
#endif
    }

    if (less(pivot,data[i]))
      swapobj(data[i],data[nData-1]);
    else {
      i++;
      //assert(i == nData-1 || less(pivot,data[i]));
      //std::cerr << "swap " << i << ":" << data[i] << " and " << (nData-1) << ":" << data[nData-1] << std::endl;
      swapobj(data[i],data[nData-1]);
    }

    //recursive calls
    const size_t n1 = i;
    const size_t n2 = nData-i-1;
    assert(n2 < nData);

    //std::cerr << "first recurse " << n1 << std::endl;
    if (n1 > 1) { //no need to sort a single entry
      quick_sort<T,Swap,Less,Equal,fallback_algo>(data,n1);
      assert(is_sorted(data,n1));
    }
    //std::cerr << "second recurse " << n2 << std::endl;
    if (n2 > 1) {
      quick_sort<T,Swap,Less,Equal,fallback_algo>(data+i+1,n2);
      assert(is_sorted(data+i+1,n2));
    }

    // std::cerr << "after recursive calls: ";
    // for (uint k=0; k < nData; k++)
    // std::cerr << data[k] << " ";
    // std::cerr << std::endl;

    assert(is_sorted(data,nData));
  }

}


template<typename T, typename Swap, typename Less, typename Equal>
void quick_sort_split_equals(T* data, const size_t nData) noexcept
{
  // std::cerr << "*** quick_sort ";
  // for (uint k=0; k < nData; k++)
  // std::cerr << k << ":" << data[k] << " ";
  // std::cerr << std::endl;

  if (nData <= QUICK_THRESH) {

    //bubble sort
    bubble_sort<T,Swap,Less>(data,nData);
  }
  else {

    const static Less less;
    const static Swap swapobj;
    using PassType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;
    PassType pivot = data[nData-1];

    size_t i=0;
    size_t j=nData-2; //not the perfect type (in view of ST)


    while(true) {

      //while((i < j) && (pivot <= data[j]))
      while((i < j) && !less(data[j],pivot))
        j--;

      //while((i < j) && (data[i] <= pivot))
      while((i < j) && !less(pivot,data[i]))
        i++;

      if (i >= j)
        break;

      assert(!less(data[i],data[j]));
      if (less(data[j],data[i])) {
        swapobj(data[i],data[j]);
#if 1
        i++;
        j--;
#endif
      }
    }


    //if (pivot < data[i])
    if (less(pivot,data[i]))
      swapobj(data[i],data[nData-1]);
    else {
      i++;
      if (less(data[nData-1],data[i]))
        swapobj(data[i],data[nData-1]);
    }

    //recursive calls
    const size_t n1 = i;
    const size_t n2 = nData-i-1;
    assert(n2 < nData);

    if (n1 > 1) { //no need to sort a single entry
      quick_sort_split_equals<T,Swap,Less,Equal>(data,n1);
      if (!is_sorted<T,Less>(data,n1))
        assert(false);
    }
    if (n2 > 1) {
      quick_sort_split_equals<T,Swap,Less,Equal>(data+i+1,n2);
      if (!is_sorted<T,Less>(data+i+1,n2))
        assert(false);
    }

    //assert(is_sorted<T,Less>(data,nData));
  }

}

template<typename T, typename Less>
inline void merge_lists(T* data, T* aux_data, size_t half, size_t nData)
{

  using RefType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;

  const static Less less;
  const size_t nData2 = nData-half;

  //now merge the lists
  const T* data1 = aux_data;
  const T* data2 = aux_data+half;

  size_t k1=0;
  size_t k2=0;
  size_t k=0;

  //while(k1 < half && k2 < nData2) {
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
}

template<typename T, typename Swap, typename Less, typename Equal>
void merge_sort(T* data, const size_t nData) noexcept
{
  //merge sort is meant for large arrays => check this for speed
  if (is_sorted<T,Less>(data, nData))
    return;

  const static Less less;

  if (nData <= MERGE_THRESH) {
    quick_sort<T,Swap,Less,Equal>(data, nData);
  }
  else {

    T* aux_data = new T[nData];
    Makros::unified_move_assign(aux_data, data, nData);
    const size_t half = nData / 2;
    const size_t nData2 = nData-half;

    merge_sort<T,Swap,Less,Equal>(aux_data,half);
    merge_sort<T,Swap,Less,Equal>(aux_data+half,nData2);

    merge_lists<T,Less>(data, aux_data, half, nData);

    // using RefType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;

    // //now merge the lists
    // const T* data1 = aux_data;
    // const T* data2 = aux_data+half;

    // size_t k1=0;
    // size_t k2=0;
    // size_t k=0;

    // //while(k1 < half && k2 < nData2) {
    // while (true) {

    // RefType d1 = data1[k1];
    // RefType d2 = data2[k2];

    // if (!less(d2,d1)) {
    // data[k] = std::move(d1);
    // k1++;
    // if (k1 >= half)
    // break;
    // }
    // else {
    // data[k] = std::move(d2);
    // k2++;
    // if (k2 >= nData2)
    // break;
    // }

    // k++;
    // }

    // //copy the rest
    // k++;
    // Makros::unified_move_assign(data+k, data1+k1, half-k1);
    // Makros::unified_move_assign(data+k, data2+k2, nData2-k2);

    delete[] aux_data;
  }
}

//recursive, won't be inline
template<typename T, typename Swap, typename Less, typename Equal>
void qmc_sort(T* data, const size_t nData) noexcept
{

  if (nData <= QUICK_THRESH) {

    //bubble sort
    bubble_sort<T,Swap,Less>(data,nData);
  }
  else {

    const static Less less;
    const static Swap swapobj;
    using PassType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;
    PassType pivot = data[nData-1];

    size_t i=0;
    size_t j=nData-2; //not the perfect type (in view of ST)

    while(true) {

      //while((i < j) && (pivot < data[j]))
      while((i < j) && less(pivot,data[j]))
        j--;

      //while((i < j) && (data[i] <= pivot))
      while((i < j) && !less(pivot,data[i]))
        i++;

      if (i >= j)
        break;

      assert(less(data[j],data[i]));
      swapobj(data[i],data[j]);
#if 0
      i++;
      j--;
#endif
    }

    if (less(pivot,data[i]))
      swapobj(data[i],data[nData-1]);
    else {
      i++;
      //assert(i == nData-1 || less(pivot,data[i]));
      //std::cerr << "swap " << i << ":" << data[i] << " and " << (nData-1) << ":" << data[nData-1] << std::endl;
      swapobj(data[i],data[nData-1]);
    }

    //recursive calls
    const size_t n1 = i;
    const size_t n2 = nData-i-1;
    assert(n2 < nData);

    //std::cerr << "n1: " << n1 << ", n2: " << n2 << std::endl;

    if (n1 > 0.667 * nData) {

      T* aux_data = new T[n1];

      Makros::unified_move_assign(aux_data, data, n1);
      const size_t half = n1 / 2;
      const size_t nData2 = n1-half;

      qmc_sort<T,Swap,Less,Equal>(aux_data,half);
      qmc_sort<T,Swap,Less,Equal>(aux_data+half,nData2);

      merge_lists<T,Less>(data, aux_data, half, n1);

      delete[] aux_data;

      if (n2 > 1) {
        qmc_sort<T,Swap,Less,Equal>(data+i+1,n2);
        //std::cerr << "AAA back from recurse" << std::endl;
        assert(is_sorted(data+i+1,n2));
      }
    }
    else if (n2 > 0.667 * nData) {

      if (n1 > 1) { //no need to sort a single entry
        qmc_sort<T,Swap,Less,Equal>(data,n1);
        //std::cerr << "BBB back from recurse" << std::endl;
        assert(is_sorted(data,n1));
      }

      T* aux_data = new T[n2];

      Makros::unified_move_assign(aux_data, data+i+1, n2);
      //std::cerr << "after move" << std::endl;
      const size_t half = n2 / 2;
      const size_t nData2 = n2-half;

      qmc_sort<T,Swap,Less,Equal>(aux_data,half);
      //std::cerr << "after first recurse" << std::endl;
      qmc_sort<T,Swap,Less,Equal>(aux_data+half,nData2);
      //std::cerr << "after second recurse" << std::endl;

      merge_lists<T,Less>(data+i+1, aux_data, half, n2);

      delete[] aux_data;
    }
    else {

      //std::cerr << "first recurse " << n1 << std::endl;
      if (n1 > 1) { //no need to sort a single entry
        qmc_sort<T,Swap,Less,Equal>(data,n1);
        assert(is_sorted(data,n1));
      }
      //std::cerr << "second recurse " << n2 << std::endl;
      if (n2 > 1) {
        qmc_sort<T,Swap,Less,Equal>(data+i+1,n2);
        assert(is_sorted(data+i+1,n2));
      }
    }

    //std::cerr << "leaving" << std::endl;
    // std::cerr << "after recursive calls: ";
    // for (uint k=0; k < nData; k++)
    // std::cerr << data[k] << " ";
    // std::cerr << std::endl;

    assert(is_sorted(data,nData));
  }
}

//interface to std::sort
template<typename T>
void std_sort(T* data, const size_t nData) noexcept
{
  std::sort(data, data+nData);
}

/**** key-value ****/

template<typename TK, typename TV, typename ValSwap, typename Less>
inline void bubble_sort_key_value(TK* key, TV* value, const size_t nData) noexcept
{
  const static ValSwap valswapobj;
  const static Less less;

  size_t last_flip = nData;

  while (last_flip > 1) {

    size_t flip = 0;

    for (size_t j=1; j < last_flip; j++) {

      if (less(key[j],key[j-1])) {

        std::swap(key[j],key[j-1]);
        valswapobj(value[j],value[j-1]);
        flip = j;
      }
    }

    last_flip = flip;
  }
}

template<typename TK, typename TV, typename ValSwap, typename Less>
inline void batch_bubble_sort_key_value(TK* key, TV* value, const size_t nData) noexcept
{
  size_t last_flip = nData;

  const static Less less;

  while (last_flip > 1) {

    //std::cerr << "*** keys: ";
    //for (uint k = 0; k < nData; k++)
    //  std::cerr << key[k] << " ";
    //std::cerr << std::endl;
    //std::cerr << "last flip: " << last_flip << std::endl;

    size_t flip = 0;

    for (size_t l=1; l < last_flip; l++) {

      //std::cerr << "l: " << l << std::endl;

      if (less(key[l],key[l-1])) {

        size_t end = l+1;
        for (; end < last_flip; end++) {
          if (!(less(key[end],key[end-1])))
            break;
        }

        //std::cerr << "reversing from " << (l-1) << " to excl. " << end << std::endl;
#ifdef USE_ROUTINES
        Routines::nontrivial_reverse<TK,SwapOp<TK> >(key+l-1, end-l+1);
        Routines::nontrivial_reverse<TV,ValSwap>(value+l-1, end-l+1);
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

template<typename T1, typename T2, typename ValSwap, typename Less>
inline void bubble_sort_key_value_prepro(T1* key, T2* value, const size_t nData) noexcept
{
  const static ValSwap valswapobj;
  const static Less less;

  const size_t thresh = 5;

  if (nData > thresh) {
    //try to ameliorate some bad cases
    for (uint i = nData-1; i >= thresh-1; i--) {
      if (less(key[i],key[0])) {
        std::swap(key[i],key[0]);
        valswapobj(value[i],value[0]);
      }
    }
  }
}

template<typename TK, typename TV, typename Less>
inline void shift_bubble_sort_key_value(TK* key, TV* value, const size_t nData) noexcept
{
  const static Less less;

  for (size_t j=1; j < nData; j++) {
    const TK curkey = key[j];
    size_t inspos = j;
    while (inspos > 0 && less(curkey,key[inspos-1]))
      inspos--;

    if (inspos != j) {
      const TV curval = std::move(value[j]);
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

template<typename TK, typename TV, typename Less> //upshift can use c++ moves, no swap involved
inline void shift_bubble_sort_key_value_keymoves(TK* key, TV* value, const size_t nData) noexcept
{
  const static Less less;

  for (size_t j=1; j < nData; j++) {
    const TK& curkey = key[j];
    size_t inspos = j;
    while (inspos > 0 && less(curkey,key[inspos-1]))
      inspos--;

    if (inspos != j) {
      const TK temp = std::move(key[j]);
      const TV curval = std::move(value[j]);
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


template<typename TK, typename TV, typename Less, typename Equal, bool stable>
inline void insertion_sort_key_value(TK* key, TV* value, const size_t nData) noexcept
{
  if (nData <= 1)
    return;

  const static Less less;
  const static Equal equal;

  if (less(key[1],key[0])) {
    std::swap(key[0],key[1]);
    std::swap(value[0],value[1]);
  }

  for (size_t i = 2; i < nData; i++) {

    size_t k = Routines::binsearch_insertpos<TK,Less,Equal>(key, key[i], i);
    if (stable) {
      while (k < i && equal(key[k],key[i]))
        k++;
    }

    if (k < i) {
      TK temp = std::move(key[i]);
      Routines::upshift_array(key, k, i, 1);
      key[k] = std::move(temp);
      TV tmp = std::move(value[i]);
      Routines::upshift_array(value, k, i, 1);
      value[k] = std::move(tmp);
    }
  }
}

//not stable and not in-place
template<typename TK, typename TV, typename Less>
inline void selection_sort_key_value(TK* key, TV* value, const size_t nData) noexcept
{
  if (nData < 2)
	return; //nothing to do

  FlexibleHeapPriorityQueueKeyValue<TK,TV,Less> heap_queue(nData);
  
  for (size_t idx = 0; idx < nData; idx++)
	heap_queue.insert(std::move(key[idx]), std::move(value[idx]));

  for (size_t idx = 0; idx < nData; idx++)
	heap_queue.extract_min(key[idx], value[idx]);
}

//in-place, but not stable
template<typename TK, typename TV, typename Less>
inline void heap_sort_key_value(TK* key, TV* value, const size_t nData) noexcept
{
  if (nData < 2)
	return; //nothing to do
	
  TK* heap = key-1;
  TV* aux_value = value-1;

  const static Less less;
  
  if (less(heap[1],heap[2])) {
    std::swap(heap[1],heap[2]);
	std::swap(aux_value[1],aux_value[2]);
  }

  //forward - insert nodes into the heap
  for(size_t idx = 3; idx <= nData; idx++) {
    //insert heap[idx] into the heap. it is already at the last position
	TK key = std::move(heap[idx]);
	const TK& ckey = key;
	
	size_t i = idx;
	
    while (i > 1) {

      const size_t parent = i >> 1; //i / 2;
      if (less(heap[parent],ckey)) {
        heap[i] = std::move(heap[parent]);
		aux_value[i] = std::move(aux_value[parent]);
	    i = parent;
      }
      else
        break;
    } 	
  
    heap[i] = std::move(key);  	  
  }

  //backward - dequeue nodes from heap
  for (size_t idx = nData; idx > 2; idx--)
  {	
	TK key = std::move(heap[idx]);
    const TK& ckey = key; 	
	TV val = std::move(aux_value[idx]);
	
	heap[idx] = std::move(heap[1]);
	
	size_t i = 1;
	
    while (true) {

      const size_t j1 = i << 1; //2*i;

      if (j1 >= idx)
        break;

      size_t better_j = j1;

	  if (j1 == idx-1) {
	    //only one child exists
	  }
      else {
        //both children exist

        const size_t j2 = j1 | 1; // j1+1;

        if (less(heap[j1],heap[j2]))
          better_j = j2;
	  }

      if (less(ckey,heap[better_j])) {
        heap[i] = std::move(heap[better_j]); 
		aux_value[i] = std::move(aux_value[better_j]);
	    i = better_j;
      }
      else
        break;
	}
	
    heap[i] = std::move(key);
	aux_value[i] = std::move(val);
  }

  std::swap(heap[1],heap[2]);
  std::swap(aux_value[1],aux_value[2]);  
}

//recursive, won't be inline
template<typename TK, typename TV, typename ValSwap,typename Less, typename Equal>
void quick_sort_key_value(TK* key, TV* value, const size_t nData) noexcept
{
  if (nData <= QUICK_THRESH) {

    //bubble sort
    bubble_sort_key_value<TK,TV,ValSwap,Less>(key,value,nData);
  }
  else {

    const static Less less;
    const static ValSwap swapobj;

    using PassType = typename std::conditional<std::is_fundamental<TK>::value || std::is_pointer<TK>::value, const TK, const TK&>::type;
    PassType pivot = key[nData-1];

    size_t i=0;
    size_t j=nData-2; //not the perfect type (in view of ST)

    while(true) {

      //while((i < j) && (pivot < key[j]))
      while((i < j) && less(pivot,key[j]))
        j--;

      //while((i < j) && (key[i] <= pivot))
      while((i < j) && !less(pivot,key[i]))
        i++;

      if (i >= j)
        break;

      assert(less(key[j],key[i]));
      std::swap(key[i],key[j]);
      swapobj(value[i],value[j]);

#if 0
      i++;
      j--;
#endif
    }

    if (less(pivot,key[i])) {
      std::swap(key[i],key[nData-1]);
      swapobj(value[i],value[nData-1]);
    }
    else {
      i++;
      assert(i == nData-1 || less(pivot,key[i]));
      std::swap(key[i],key[nData-1]);
      swapobj(value[i],value[nData-1]);
    }

    //recursive calls
    const size_t n1 = i;
    const size_t n2 = nData-i-1;
    assert(n2 < nData);

    if (n1 > 1) { //no need to sort a single entry
      quick_sort_key_value<TK,TV,ValSwap,Less,Equal>(key,value,n1);
      assert(is_sorted(key,n1));
    }
    if (n2 > 1) {
      quick_sort_key_value<TK,TV,ValSwap,Less,Equal>(key+i+1,value+i+1,n2);
      assert(is_sorted(key+i+1,n2));
    }

    assert(is_sorted(key,nData));
  }
}

template<typename TK, typename TV>
inline void merge_lists(TK* key, TV* value, TK* aux_key, TV* aux_value, const size_t half, const size_t nData)
{

  using KeyRefType = typename std::conditional<std::is_fundamental<TK>::value || std::is_pointer<TK>::value, const TK, const TK&>::type;

  const size_t nData2 = nData-half;

  const TK* key1 = aux_key;
  const TK* key2 = aux_key+half;
  const TV* value1 = aux_value;
  const TV* value2 = aux_value+half;

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
}

//recursive, won't be inline
template<typename TK, typename TV, typename Swap, typename Less, typename Equal>
void merge_sort_key_value(TK* key, TV* value, const size_t nData) noexcept
{
  //merge sort is meant for large arrays => check this for speed
  if (is_sorted(key, nData))
    return;

  if (nData <= MERGE_THRESH)
    quick_sort_key_value<TK,TV,Swap,Less,Equal>(key, value, nData);
  else {

    using KeyRefType = typename std::conditional<std::is_fundamental<TK>::value || std::is_pointer<TK>::value, const TK, const TK&>::type;

    const size_t half = nData / 2;
    const size_t nData2 = nData-half;

    TK* aux_key = new TK[nData];
    TV* aux_value = new TV[nData];

    Makros::unified_move_assign(aux_key, key, nData);
    Makros::unified_move_assign(aux_value, value, nData);

    merge_sort_key_value(aux_key,aux_value,half);
    merge_sort_key_value(aux_key+half,aux_value+half,nData2);

    //now merge the lists
    merge_lists<TK,TV>(key, value, aux_key, aux_value, half, nData);


    // const TK* key1 = aux_key;
    // const TK* key2 = aux_key+half;
    // const TV* value1 = aux_value;
    // const TV* value2 = aux_value+half;

    // size_t k1=0;
    // size_t k2=0;
    // size_t k=0;

    // while (true) {

    // const KeyRefType d1 = key1[k1];
    // const KeyRefType d2 = key2[k2];
    // if (d1 <= d2) {
    // key[k] = std::move(d1);
    // value[k] = std::move(value1[k1]);
    // k1++;
    // if (k1 >= half)
    // break;
    // }
    // else {
    // key[k] = std::move(d2);
    // value[k] = std::move(value2[k2]);
    // k2++;
    // if (k2 >= nData2)
    // break;
    // }

    // k++;
    // }

    // //copy the rest
    // k++;
    // Makros::unified_move_assign(key+k, key1+k1, half-k1);
    // Makros::unified_move_assign(value+k, value1+k1, half-k1);
    // Makros::unified_move_assign(key+k, key2+k2, nData2-k2);
    // Makros::unified_move_assign(value+k, value2+k2, nData2-k2);

    delete[] aux_key;
    delete[] aux_value;
  }
}

//recursive, won't be inline
template<typename TK, typename TV, typename ValSwap, typename Less, typename Equal>
void qmc_sort_key_value(TK* key, TV* value, const size_t nData) noexcept
{

  if (nData <= QUICK_THRESH) {

    //bubble sort
    bubble_sort_key_value<TK,TV,ValSwap,Less>(key,value,nData);
  }
  else {

    const static Less less;
    const static ValSwap swapobj;
    using PassType = typename std::conditional<std::is_fundamental<TK>::value || std::is_pointer<TK>::value, const TK, const TK&>::type;
    PassType pivot = key[nData-1];

    size_t i=0;
    size_t j=nData-2; //not the perfect type (in view of ST)

    while(true) {

      while((i < j) && less(pivot,key[j]))
        j--;

      while((i < j) && !less(pivot,key[i]))
        i++;

      if (i >= j)
        break;

      assert(less(key[j],key[i]));
      std::swap(key[i],key[j]);
      swapobj(value[i],value[j]);

#if 0
      i++;
      j--;
#endif
    }

    if (less(pivot,key[i])) {
      std::swap(key[i],key[nData-1]);
      swapobj(value[i],value[nData-1]);
    }
    else {
      i++;
      //assert(i == nData-1 || less(pivot,data[i]));
      //std::cerr << "swap " << i << ":" << data[i] << " and " << (nData-1) << ":" << data[nData-1] << std::endl;
      std::swap(key[i],key[nData-1]);
      swapobj(value[i],value[nData-1]);
    }

    //recursive calls
    const size_t n1 = i;
    const size_t n2 = nData-i-1;
    assert(n2 < nData);

    //std::cerr << "n1: " << n1 << ", n2: " << n2 << std::endl;

    if (n1 > 0.667 * nData) {

      TK* aux_key = new TK[n1];
      TV* aux_value = new TV[n1];

      Makros::unified_move_assign(aux_key, key, n1);
      Makros::unified_move_assign(aux_value, value, n1);

      const size_t half = n1 / 2;
      const size_t nData2 = n1-half;

      qmc_sort_key_value<TK,TV,ValSwap,Less,Equal>(aux_key,aux_value,half);
      qmc_sort_key_value<TK,TV,ValSwap,Less,Equal>(aux_key+half,aux_value+half,nData2);

      //now merge the lists
      merge_lists<TK,TV>(key, value, aux_key, aux_value, half, n1);


      delete[] aux_key;
      delete[] aux_value;

      if (n2 > 1) {
        qmc_sort_key_value<TK,TV,ValSwap,Less,Equal>(key+i+1,value+i+1,n2);
        //std::cerr << "AAA back from recurse" << std::endl;
        assert(is_sorted(key+i+1,n2));
      }
    }
    else if (n2 > 0.667 * nData) {

      if (n1 > 1) { //no need to sort a single entry
        qmc_sort_key_value<TK,TV,ValSwap,Less,Equal>(key,value,n1);
        //std::cerr << "BBB back from recurse" << std::endl;
        assert(is_sorted(key,n1));
      }

      TK* aux_key = new TK[n2];
      TV* aux_value = new TV[n2];

      Makros::unified_move_assign(aux_key, key+i+1, n2);
      Makros::unified_move_assign(aux_value, value+i+1, n2);

      //std::cerr << "after move" << std::endl;
      const size_t half = n2 / 2;
      const size_t nData2 = n2-half;

      qmc_sort_key_value<TK,TV,ValSwap,Less,Equal>(aux_key,aux_value,half);
      //std::cerr << "after first recurse" << std::endl;
      qmc_sort_key_value<TK,TV,ValSwap,Less,Equal>(aux_key+half,aux_value+half,nData2);
      //std::cerr << "after second recurse" << std::endl;

      merge_lists<TK,TV>(key+i+1, value+i+1, aux_key, aux_value, half, n2);

      delete[] aux_key;
      delete[] aux_value;
    }
    else {

      //std::cerr << "first recurse " << n1 << std::endl;
      if (n1 > 1) { //no need to sort a single entry
        qmc_sort_key_value<TK,TV,ValSwap,Less,Equal>(key,value,n1);
        assert(is_sorted(key,n1));
      }
      //std::cerr << "second recurse " << n2 << std::endl;
      if (n2 > 1) {
        qmc_sort_key_value<TK,TV,ValSwap,Less,Equal>(key+i+1,value+i+1,n2);
        assert(is_sorted(key+i+1,n2));
      }
    }

    //std::cerr << "leaving" << std::endl;
    // std::cerr << "after recursive calls: ";
    // for (uint k=0; k < nData; k++)
    // std::cerr << data[k] << " ";
    // std::cerr << std::endl;

    assert(is_sorted(key,nData));
  }

  TODO("");
}


/**** index ****/

template<typename T, typename ST>
inline void aux_index_bubble_sort(const T* data, ST* indices, const ST nData) noexcept
{
  size_t last_flip = nData;

  while (last_flip > 1) {

    size_t flip = 0;

    for (size_t l=1; l < last_flip; l++) {

      if (data[indices[l]] < data[indices[l-1]]) {
        std::swap(indices[l],indices[l-1]);
        flip = l;
      }
    }

    last_flip = flip;
  }

}

template<typename T, typename ST>
inline void index_bubble_sort(const T* data, ST* indices, const ST nData) noexcept
{
  for (uint i=0; i < nData; i++)
    indices[i] = i;

  aux_index_bubble_sort(data, indices, nData);
}

template<typename T, typename ST>
inline void aux_index_batch_bubble_sort(const T* data, ST* indices, ST nData) noexcept
{
  const static std::less<T> less;

  size_t last_flip = nData;

  while (last_flip > 1) {

    size_t flip = 0;

    for (size_t l=1; l < last_flip; l++) {

      if (less(data[indices[l]],data[indices[l-1]])) {

        //std::cerr << "found pair at " << l << std::endl;

        size_t end = l+1;
        for (; end < last_flip; end++) {

          if (!less(data[indices[end]],data[indices[end-1]]))
            break;
        }

#ifdef USE_ROUTINES
        Routines::nontrivial_reverse<T,SwapOp<ST> >(indices+l-1, end-l+1);
#else
        std::reverse(indices+l-1, indices+end);
#endif

        flip = end-1;
        l = end-2; //the last swapped index needs to be compared again (l is incremented by the loop)
      }
    }

    last_flip = flip;
  }
}

template<typename T, typename ST>
inline void index_batch_bubble_sort(const T* data, ST* indices, ST nData) noexcept
{
  for (ST i=0; i < nData; i++)
    indices[i] = i;

  aux_index_batch_bubble_sort(data,indices,nData);
}

template<typename T, typename ST, bool stable>
inline void aux_index_insertion_sort(T* data, ST* indices, const ST nData) noexcept
{
  if (nData <= 1)
    return;

  const static std::less<T> less;
  const static std::equal_to<T> equal;

  if (less(data[indices[1]],data[indices[0]]))
    std::swap(indices[0],indices[1]);

  for (ST i = 2; i < nData; i++) {

    ST k = Routines::index_binsearch_insertpos<T,ST>(data, data[i], indices, i);
    //std::cerr << "key: " << data[i] << ", found pos: " << k << std::endl;
    if (stable) {
      while (k < i && equal(data[indices[k]],data[indices[i]]))
        k++;
    }

    if (k < i) {
      ST temp = std::move(indices[i]);
      Routines::upshift_array(indices, k, i, 1);
      indices[k] = std::move(temp);
    }
  }
}

template<typename T, typename ST, bool stable>
inline void index_insertion_sort(T* data, ST* indices, const ST nData) noexcept
{
  for (ST i=0; i < nData; i++)
    indices[i] = i;

  aux_index_insertion_sort<T,ST,stable>(data,indices,nData);
}

template<typename T, typename ST>
inline void index_selection_sort(const T* data, ST* indices, ST nData) noexcept
{
  for (ST i=0; i < nData; i++)
    indices[i] = i;

  if (nData < 2)
    return;	  

  FlexibleHeapPriorityQueueIndexed<T,ST> heap_queue(nData,data);

  for (size_t idx = 0; idx < nData; idx++)
    heap_queue.insert(std::move(indices[idx]));	  
 
  //std::cerr << "proceed" << std::endl;
  for (size_t idx = 0; idx < nData; idx++)
    heap_queue.extract_min(indices[idx]);
}

template<typename T, typename ST>
inline void index_heap_sort(const T* data, ST* indices, ST nData) noexcept
{
  for (ST i=0; i < nData; i++)
    indices[i] = i;

  if (nData < 2)
    return;	  

  ST* heap = indices-1;

  if (data[heap[1]] < data[heap[2]])
    std::swap(heap[1],heap[2]);

  //forward - insert nodes into the heap
  for(size_t idx = 3; idx <= nData; idx++) {
    //insert heap[idx] into the heap. it is already at the last position
    const ST cidx = indices[idx];
	const T curkey = data[cidx];
	
	size_t i = idx;
	
    while (i > 1) {

      const size_t parent = i >> 1; //i / 2;
      if (data[heap[parent]] < curkey) {
        heap[i] = heap[parent];
	    i = parent;
      }
      else
        break;
    } 	
  
    heap[i] = cidx;
  }	  

  //backward - dequeue nodes from heap
  for (size_t idx = nData; idx > 2; idx--)
  {	
	const ST key = std::move(heap[idx]);
	heap[idx] = std::move(heap[1]);
	const T ckey = data[key];
	
	size_t i = 1;
	
    while (true) {

      const size_t j1 = i << 1; //2*i;

      if (j1 >= idx)
        break;

      size_t better_j = j1;

	  if (j1 == idx-1) {
	    //only one child exists
	  }
      else {
        //both children exist

        const size_t j2 = j1 | 1; // j1+1;

        if (data[heap[j1]] < data[heap[j2]])
          better_j = j2;
	  }

      if (ckey < data[heap[better_j]]) {
        heap[i] = std::move(heap[better_j]); 
	    i = better_j;
      }
      else
        break;
	}
	
    heap[i] = key;
  }

  std::swap(heap[1],heap[2]);  
}


//recursive, won't be inline
template<typename T, typename ST>
void aux_index_quick_sort(const T* data, ST* indices, const ST nData) noexcept
{

  //std::cerr << "nData: " << nData << std::endl;

  if (nData <= QUICK_THRESH) {

    aux_index_bubble_sort<T,ST>(data,indices,nData);
  }
  else {

    using PassType = typename std::conditional<std::is_fundamental<T>::value || std::is_pointer<T>::value, const T, const T&>::type;
    PassType pivot = data[indices[nData-1]];

    ST i=0;
    ST j=nData-2; //not the perfect type (in view of ST)

    while (true) {

      while(i < j && data[indices[i]] <= pivot)
        i++;

      while(j > i && data[indices[j]] > pivot)
        j--;

      if (i >= j)
        break;

      std::swap(indices[i],indices[j]);
#if 1
      i++;
      j--;
#endif
    }

    if (pivot < data[indices[i]])
      std::swap(indices[i],indices[nData-1]);
    else {
      i++;
      std::swap(indices[i],indices[nData-1]);
    }

    //recursive calls
    const ST n1 = i;
    const ST n2 = nData-i-1;

    if (n1 > 1) { //no need to sort a single entry
      aux_index_quick_sort<T,ST>(data,indices,n1);
    }
    if (n2 > 1) {
      aux_index_quick_sort<T,ST>(data,indices+i+1,n2);
    }
  }

}

template<typename T, typename ST>
void index_quick_sort(const T* data, ST* indices, const ST nData) noexcept
{
  //std::cerr << "sorting " << nData << " entries" << std::endl;

  for (uint i=0; i < nData; i++)
    indices[i] = i;

  aux_index_quick_sort(data,indices,nData);
}

//recursive, won't be inline
template<typename T, typename ST>
void aux_index_merge_sort(const T* data, ST* indices, const ST nData) noexcept
{
  //merge sort is meant for large arrays => check this for speed
  if (is_sorted(data, nData))
    return;

  if (nData <= MERGE_THRESH) {
    aux_index_quick_sort(data,indices,nData);
  }
  else {

    const ST half = nData / 2;
    const ST nData2 = nData-half;

    ST* aux_indices = new ST[nData];

    Makros::unified_move_assign(aux_indices, indices, nData);

    aux_index_merge_sort(data,aux_indices,half);
    aux_index_merge_sort(data,aux_indices+half,nData2);

    const ST* index1 = aux_indices;
    const ST* index2 = aux_indices+half;

    ST k1=0;
    ST k2=0;
    ST k=0;

    //while(k1 < half && k2 < nData2) {
    while (true) {

      const ST idx1 = index1[k1];
      const ST idx2 = index2[k2];
      if (data[idx1] <= data[idx2]) {
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

//this is still fixed to memcpy
template<typename T, typename ST>
void aux_index_merge_sort_4split(const T* data, ST* indices, const ST nData) noexcept
{
  //std::cerr << "nData: "<< nData << std::endl;

  if (nData <= 8 /*8*/) {

    //bubble sort
    for (uint k=0; k < nData-1; k++) {

      //std::cerr << "k: " << k << std::endl;

      for (uint l=0; l < nData-1-k; l++) {
        //std::cerr << "l: " << l << std::endl;

        //std::cerr << "index1: " << indices[l] << std::endl;
        //std::cerr << "index2: " << indices[l+1] << std::endl;

        if (data[indices[l]] > data[indices[l+1]])
          std::swap(indices[l],indices[l+1]);
      }
    }
  }
  else if (nData <= 32) {
    aux_index_merge_sort(data,indices,nData);
  }
  else {

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

    aux_index_merge_sort(data,index[0].direct_access(),nData1to3);

    index[1].resize(nData1to3);
    memcpy(index[1].direct_access(),indices+nData1to3,nData1to3*sizeof(ST));

    aux_index_merge_sort(data,index[1].direct_access(),nData1to3);

    index[2].resize(nData1to3);
    memcpy(index[2].direct_access(),indices+2*nData1to3,nData1to3*sizeof(ST));

    aux_index_merge_sort(data,index[2].direct_access(),nData1to3);

    index[3].resize(nData4);
    memcpy(index[3].direct_access(),indices+3*nData1to3,nData4*sizeof(ST));

    aux_index_merge_sort(data,index[3].direct_access(),nData4);

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

          if (arg_best == MAX_UINT || cur_val < best) {
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

template<typename T, typename ST>
void index_merge_sort(const T* data, ST* indices, const ST nData) noexcept
{
  //std::cerr << "sorting " << nData << " entries" << std::endl;

  for (ST i=0; i < nData; i++)
    indices[i] = i;

  aux_index_merge_sort(data,indices,nData);
}





#endif
