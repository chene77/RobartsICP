//
// Elvis Chen
// chene@cs.queensu.ca
//
// Department of Computing and Information Science
// Queen's University, Kingston, Ontario, Canada
//
// Feb. 20, 2000
//

//
// Filename:  qsort.h
//

// Initial implementation of QuickSort
//
// not that we are performing an indexing sort; the array
// itself is not sorted.
//

#ifndef __QSORT_H__
#define __QSORT_H__

#include "matrix.h"

template<class T>
inline void SWAP( T &a, T &b ) 
{
  T temp = a;
  a = b;
  b = temp;
}

template< class T >
Vec<Subscript> qsort( const Vec<T> &arr ) 
{
  Subscript M = 7;
  Subscript NSTACK = 50;
  Subscript n = arr.dim();
  
  Subscript i, indxt, ir = n - 1, j, k, l = 0, jstack = 0;
  T a;
  
  Vec<Subscript> indx( n );
  Vec<Subscript> istack( NSTACK );
  
  for (j = 0; j < n; j++) indx[j] = j;
    
  for (;;) {
    if (ir-l < M) {
      // insertion sort when subarray small enough.
      for (j = l+1; j <= ir; j++) {
	indxt = indx[j];
	a = arr[indxt];
	for (i = j - 1; i >= l; i--) {
	  if (arr[indx[i]] <= a) break;
	  indx[i+1] = indx[i];
	}
	indx[i+1] = indxt;
      }
      if (jstack == 0) break;
      ir = istack[jstack--];
      l = istack[jstack--];
    }
    else {
      k=(l+ir) >> 1;
      SWAP(indx[k], indx[l+1]);
      if (arr[indx[l]] > arr[indx[ir]]) {
	SWAP(indx[l], indx[ir]);
      }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	SWAP(indx[l+1], indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
	SWAP(indx[l], indx[l+1]);
      }
      i = l+1;
      j = ir;
      indxt = indx[l+1];
      
      
      a = arr[indxt];
      for (;;){
	do i++;
	while (arr[indx[i]] < a);
	
	do j--;
	while (arr[indx[j]] > a);
	
	if (j < i) break;
	SWAP(indx[i], indx[j]);
      }
      
      indx[l+1] = indx[j];
      indx[j] = indxt;
      jstack += 2;
      
      assert( jstack <= NSTACK );
      
      /*
      if (jstack > NSTACK) 
	Elerror("NSTACK too small in qsort()");
	*/
      
      if (ir - i + 1 >= j-l) {
	istack[jstack] = ir;
	istack[jstack-1] = i;
	ir = j-1;
      }
      
      else {
	istack[jstack] = j-1;
	istack[jstack-1] = l;
	l = i;
      }
    }
  }
      
  return indx;
  
}


#endif // of __QSORT_H__
