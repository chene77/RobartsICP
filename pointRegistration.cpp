//
// Elvis Chen
// chene@cs.queensu.ca
//
// Department of Computing and Information Science
// Queen's University, Kingston, Ontario, Canada
//
// Feb. 16, 2000
//

//
// Filename:  matrix.h
//

// Initial implementation of a matrix class, based on TNT (Template
// Numerical Toolkit, http://math.nist.gov/tnt).
//
//
// C compatible matrix:  row-oriented, 0-based[i][j] indexing
//
// Numerical Recipes compatible matrix:  row-oriented, 1-based(i,j) indexing
//

// This is actually a subset/change of TNT, since we don't really
// care about Fortran-style (column-oriented, 1-based) indexing.  
// However, we would like to have ROW-oriented, 1-based indexing
// so we can implement Numerical Recipes with ease.


#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <cassert>
#include <iostream>
#include <strstream>
#include <cmath>
#include <cstdlib>

// Quote from the web page:
//
// A good deal of the segmentation faults produced by numerical
// codes are related to indexing a vector or matrix past its 
// pre-allocated bounds. By default, TNT verifies that the index 
// used to access vector elements is always within the vector's bounds.
// (Just like Java does.) This can be turned off by using the
// following compile flag 
//
//                         -DNO_BOUNDS_CHECK
//
// which removes all performance penalties incurred for performing
// this check. (Typically done only after final debugging.) 

#define BOUNDS_CHECK
#ifdef NO_BOUNDS_CHECK
#undef BOUNDS_CHECK
#endif


//---------------------------------------------------------------------
// Define the data type used for matrix and vector Subscripts.
// This will default to "int", but it can be overriden at compile time,
// e.g.
//
//      g++ -DMATRIX_SUBSCRIPT_TYPE='unsinged long' ...
//
//---------------------------------------------------------------------

#ifndef MATRIX_SUBSCRIPT_TYPE
#define MATRIX_SUBSCRIPT_TYPE int
#endif

typedef MATRIX_SUBSCRIPT_TYPE Subscript;


template <class T> class Vec 
{
public:
  typedef Subscript size_type;
  typedef T         value_type;
  typedef T         element_type;
  typedef T*        pointer;
  typedef T*        iterator;
  typedef T&        reference;
  typedef const T*  const_iterator;
  typedef const T&  const_reference;
  
protected:
  T *_v;  // c-style indexing, 0-based
  T *_v1; // row-major, 1-based
  Subscript _size;
  
  // internal helper function to create the array 
  // of row pointers
  
  void initialize( Subscript N ) 
  {
    assert( _v == NULL );
  
    _v = new T[N];
    assert( _v != NULL );
    
    _v1 = _v - 1;

    _size = N;
  }
  
  void copy( const T *v )
  {
    for (Subscript i = 0; i < _size; i++) _v[i] = v[i];
  }
  
  void set( const T &val )
  {    
    for (Subscript i = 0; i < _size; i++) _v[i] = val;
  }
  
  void destroy() 
  {
    if ( _v == NULL) return;
    
    delete [] ( _v );
    _v = NULL;
    _v1 = NULL;
  }
  
public:
  // access
  iterator begin() { return _v; }
  iterator end() { return (_v + _size); }
  const iterator begin() const { return _v; }
  const iterator end() const{ return ( _v + _size ); }
  
  // destructor
  ~Vec() { destroy(); }
  
  // constructors
  Vec() : _v(0), _v1(0), _size(0) {};
  
  Vec( const Vec<T> &A ) : _v(0), _v1(0), _size(0) 
  {
    initialize( A._size );
    copy( A._v );
  };
    
  Vec( Subscript N, const T &value = T(0) ) : _v(0), _v1(0), _size(0) 
  {
    initialize( N );
    set( value );
  };
    
  Vec( Subscript N, const T *v ): _v(0), _v1(0),_size(0) 
  {
    initialize( N );
    copy( v );
  };
    
  Vec( Subscript N, char *s ) : _v(0), _size(0) 
  {
    initialize( N );
    std::istrstream ins(s);
    
    Subscript i;
    
    for (i = 0; i < N; i++) ins >> _v[i];
  };
  
  
  //
  // methods
  //
  Vec<T> &newsize( Subscript N )
  {
    if ( _size == N ) return *this;
    
    destroy();
    initialize( N );
    
    return *this;
  };
  
  //
  // assignments
  // 
  Vec<T>& operator= (const Vec<T> &A) 
  {
    if ( _v == A._v ) return *this;
    
    if ( _size == A._size ) // no need to re-allocate
      copy( A._v );
    else {
      destroy();
      initialize( A._size );
      copy( A._v );
    }
    
    return *this;
  };
    
  Vec<T>& operator= (const T &scalar) 
  {
    set(scalar);
    return *this;
  };
  
  Vec<T>& operator+=( const Vec<T> &A ) 
  {
    assert( _size == A._size );
    
    for (Subscript i = 0; i < _size; i++) _v[i] += A[i];
    return *this;
  };
  
  template<class S>
  Vec<T>& operator+=( const S &val ) 
  {
    for (Subscript i = 0; i < _size; i++) _v[i] += (T)val;
    return *this;
  };
  
  Vec<T>& operator-=( const Vec<T> &A ) 
  {
    assert( _size == A._size );
    
    for (Subscript i = 0; i < _size; i++) _v[i] -= A[i];
    return *this;
  };

  template<class S>
  Vec<T>& operator-=( const S &val ) 
  {
    for (Subscript i = 0; i < _size; i++) _v[i] -= (T)val;
    return *this;
  };

  Vec<T>& operator*=( const Vec<T> &A ) 
  {
    assert( _size == A._size );
    
    for (Subscript i = 0; i < _size; i++) _v[i] *= A[i];
    return *this;
  };

  template<class S>
  Vec<T>& operator*=( const S &val ) 
  {
    for (Subscript i = 0; i < _size; i++) _v[i] *= (T)val;
    return *this;
  };

  Vec<T>& operator/=( const Vec<T> &A ) 
  {
    assert( _size == A._size );
    
    for (Subscript i = 0; i < _size; i++) _v[i] /= A[i];
    return *this;
  };

  template<class S>
  Vec<T>& operator/=( const S &val ) 
  {
    for (Subscript i = 0; i < _size; i++) _v[i] /= (T)val;
    return *this;
  };
  
  inline Subscript dim() const 
  {
    return _size;
  };
    
  inline Subscript size() const 
  {
    return _size;
  };
    
  inline reference operator() ( Subscript i )
  {
#ifdef BOUND_CHECK
    assert( 1 <= i );
    assert( i <= _size );
#endif
  
    return _v1[ i ];
  };
    
  inline const_reference operator() ( Subscript i ) const 
  {
#ifdef BOUND_CHECK
    assert( 1 <= i );
    assert( i <= _size );
#endif
    
    return _v1[ i ];
  };
    
  inline reference operator[] ( Subscript i )
  {
#ifdef BOUNDS_CHECK
    assert( 0 <= i );
    assert( i < _size );
#endif

    return _v[i];
  };
    
  inline const_reference operator[] ( Subscript i ) const 
  {
#ifdef BOUNDS_CHECK
    assert( 0 <= i );
    assert( i < _size );
#endif

    return _v[i];
  };
  
  // min and max values of the given vector
  inline T min() 
  {
    T tmp = _v[0];
    
    for (Subscript i = 1; i < _size; i++) if ( _v[i] < tmp ) tmp = _v[i];
    return (tmp);
  };
  
  inline T max() 
  {
    T tmp = _v[0];
    
    for (Subscript i = 1; i < _size; i++) if ( _v[i] > tmp ) tmp = _v[i];
    return (tmp) ;
  };
};



template <class T> class Matrix 
{
public:
  typedef Subscript size_type;
  typedef T         value_type;
  typedef T         element_type;
  typedef T*        pointer;
  typedef T*        iterator;
  typedef T&        reference;
  typedef const T*  const_iterator;
  typedef const T&  const_reference;
  

protected:
  Subscript size_x;
  Subscript size_y;
  Subscript totalSize; // total size
  T *_v;
  T **_row;  // c-style, 0-based

  T *_v1;
  T **_row1; // row-major, 1-based
  
  // internal helper function to create the array 
  // of row pointers

  void initialize( Subscript X, Subscript Y ) 
  {
    totalSize = X * Y;
    size_x = X;
    size_y = Y;
    
    _v = new T[ totalSize ];
    _row = new T*[ X ];
    _row1 = new T*[ X ];
    
    assert( _v != NULL );
    assert( _row != NULL );
    assert( _row1 != NULL );
    
    T *p = _v;
    _v1 = _v - 1;
    
    for (Subscript i = 0; i < X; i++) {
      _row[i] = p;
      _row1[i] = p - 1;
      
      p += Y;
    }

    _row1 --;  // compensate for 1-based offset
  };
  
  // copy two matrices
  void copy( const T *v ) 
  {
    for (Subscript i = 0; i < (size_x * size_y); i++) _v[i] = v[i];
  };
    
  // Initial the matrix with a given value
  void set( const T& val )
  {
    for (Subscript i = 0; i < (size_x * size_y); i++) _v[i] = val;
  };
    
  // helper function for destructor
  void destroy() 
  {
    // do nothing, if no memory has been previously allocated
    if ( _v == NULL ) return;
    
    // de-allocate the memory
    if ( _v != NULL ) delete [] ( _v );
    if ( _row != NULL ) delete [] ( _row );
        
    // return _row1 back to the original value
    _row1 ++;
    if ( _row1 != NULL) delete [] ( _row1 );
    
    _v = NULL;
    _row = NULL;
    _row1 = NULL;
    
  };
  
public:

  operator T**() { return _row; }
  operator T**() const { return _row; }

  Subscript size() const { return totalSize; }


  // constructors
  Matrix() : size_x(0), size_y(0), totalSize(0), 
    _v(0), _row(0), _v1(0), _row1(0) {};
  
  Matrix( const Matrix<T> &A ) 
  {
    initialize( A.size_x, A.size_y );
    copy( A._v );
  };
    
  Matrix( Subscript X, Subscript Y, const T& value = T(0) ) 
  {
    initialize( X, Y );
    set( value );
  };
    
  Matrix( Subscript X, Subscript Y, const T *v ) 
  {
    initialize( X, Y );
    copy( v );
  };
    
  Matrix(Subscript X, Subscript Y, char *s )
  {
    initialize( X, Y );
    std::istrstream ins(s);
    
    Subscript i, j;
    
    for (i = 0; i < X; i++)
      for (j = 0; j < Y; j++)
	ins >> _row[i][j];
  };
    

  //
  // descructor
  //
  ~Matrix() 
  {
    destroy();
  };
  
  // reallocating
  Matrix<T>& newsize( Subscript X, Subscript Y )
  {
    if ( num_rows() == X && num_cols() == Y ) return *this;
    
    destroy();
    initialize( X, Y );
    
    return *this;
  };
    
  Matrix<T>& newsize( Subscript X, Subscript Y, const T &value )
  {
    if ( num_rows() == X && num_cols() == Y ) {
      set( value );
      return *this;
    }
    
    destroy();
    initialize( X, Y );
    set( value );
    
    return *this;
  };
  
  //
  // assignments
  //
  Matrix<T>& operator= (const Matrix<T> &A) 
  {
    if ( _v == A._v ) return *this;
    
    if ( size_x == A.size_x && size_y == A.size_y ) // no need to re-allocate
      copy( A._v );
    else {
      destroy();
      initialize( A.size_x, A.size_y );
      copy( A._v );
    }
    
    return *this;
  };
    
  Matrix<T>& operator= (const T& scalar) 
  {
    set(scalar);
    return *this;
  };
    
  Matrix<T>& operator+= (const Matrix<T> &A) 
  {
    assert ( (size_x == A.size_x) && (size_y == A.size_y) );
    for (Subscript i = 0; i < totalSize; i++) _v[i] += A._v[i];
    return *this;
  };
  
  template<class S>
  Matrix<T>& operator+= (const S &scalar) 
  {
    T val = (T)scalar;
    for (Subscript i = 0; i < totalSize; i++) _v[i] += val;
    return *this;
  };
    
  Matrix<T>& operator-= (const Matrix<T> &A) 
  {
    assert ( (size_x == A.size_x) && (size_y == A.size_y) );
    for (Subscript i = 0; i < totalSize; i++) _v[i] -= A._v[i];
    return *this;
  };

  template<class S>
  Matrix<T>& operator-= (const S &scalar) 
  {
    T val = (T)scalar;
    for (Subscript i = 0; i < totalSize; i++) _v[i] -= val;
    return *this;
  };
  
  template<class S>
  Matrix<T>& operator*= (const S &scalar) 
  { 
    T val = (T)scalar;
    for (Subscript i = 0; i < totalSize; i++) _v[i] *= val;
    return *this;
  };
    
  template<class S>
  Matrix<T>& operator/= (const S &scalar) 
  {
    T val = (T)scalar;
    for (Subscript i = 0; i < totalSize; i++) _v[i] /= val;
    return *this;
  };
  
  Subscript dim(Subscript d) const 
  {
#ifdef BOUNDS_CHECK
    assert( d >= 1 );
    assert( d <= 2 );
#endif
    return (d == 1) ? size_x : ((d==2) ? size_y : 0);
  };
    
  Subscript num_rows() const { return size_x; }
  Subscript num_cols() const { return size_y; }
  
  inline T* operator[] (Subscript i) 
  {
#ifdef BOUNDS_CHECK
    assert( 0 <= i );
    assert( i < size_x );
#endif
    return _row[i];
  };
    
  inline const T* operator[] (Subscript i) const 
  {
#ifdef BOUNDS_CHECK
    assert( 0 <= i );
    assert( i < size_x );
#endif
    return _row[i];
  };
    
  inline reference operator() (Subscript i) 
  {
#ifdef BOUNDS_CHECK
    assert( 1 <= i );
    assert( i <= totalSize );
#endif
    return _v1[i];
  };
    
  inline const_reference operator() (Subscript i) const
  {
#ifdef BOUNDS_CHECK
    assert( 1 <= i );
    assert( i <= totalSize );
#endif
    return _v1[i-1];
  };
    
  inline reference operator() (Subscript i, Subscript j) 
  {
#ifdef BOUNDS_CHECK
    assert( 1 <= i );
    assert( i <= size_x );
    assert( 1 <= j);
    assert( j <= size_y );
#endif
    return _row1[i][j];
  };
  
  inline const_reference operator() (Subscript i, Subscript j) const 
  {
#ifdef BOUNDS_CHECK
    assert( 1 <= i );
    assert( i <= size_x );
    assert( 1 <= j );
    assert( j <= size_y );
#endif
    return _row1[i][j];
  };
  
  // min/max value of the given matrix
  inline T min() 
  {
    T tmp = _v[0];
    
    for (Subscript i = 1; i < totalSize; i++) if ( _v[i] < tmp ) tmp = _v[i];
    
    return (tmp);
  };
  
  inline T max() 
  {
    T tmp = _v[0];
    
    for (Subscript i = 1; i < totalSize; i++) if ( _v[i] > tmp ) tmp = _v[i];
    
    return (tmp);
  };
};

//
// I/O operation
//

//
// I/O for Vec
//
template <class T>
inline std::ostream& operator<< (std::ostream &s, const Vec<T> &A) 
{
  Subscript N = A.dim();
  
  s << N << std::endl;
  
  for (Subscript i = 0; i < N; i++) s << A[i] << " " << std::endl;
  
  s << std::endl;
  
  return s;
}

template <class T>
inline std::istream & operator>> (std::istream &s, Vec<T> &A) 
{
  Subscript N;
  
  s >> N;
  
  if ( !(N == A.size()) ){
    A.destroy();
    A.initialize(N);
  }
  
  for (Subscript i = 0; i < N; i++) s >> A[i];
  
  return s;
}

//
// I/O for Matrix
//
template <class T>
inline std::ostream& operator<< ( std::ostream &s,
				  const Matrix<T> &A ) 
{
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  s << X << " " << Y << std::endl;
  
  for (Subscript i = 0; i < X; i++) {
    for (Subscript j = 0; j < Y; j++) {
      s << A[i][j] << " ";
    }
    s << std::endl;
  }
  
  return s;
}


template <class T> 
inline std::istream& operator>> (std::istream &s,
				 Matrix<T> &A )
{
  Subscript X, Y;
  
  s >> X >> Y;
  
  if ( !(X == A.num_rows() && Y == A.num_cols()) ){
    A.newsize( X, Y );
  }
  
  
  for (Subscript i = 0; i < X; i++) 
    for (Subscript j = 0; j < Y; j++) 
      s >> A[i][j];
  
  return s;
}


//
// basic vector algorithms
//
template <class T>
inline Vec<T> operator+ (const Vec<T> &A, const Vec<T> &B) 
{
  Subscript N = A.dim();
  
  assert( N == B.dim() );
  
  Vec<T> tmp(N);
  
  for (Subscript i = 0; i < N; i++) tmp[i] = A[i] + B[i];
  
  return tmp;
}

template <class T>
inline Vec<T> operator- (const Vec<T> &A, const Vec<T> &B) 
{
  Subscript N = A.dim();
  
  assert( N == B.dim() );
  
  Vec<T> tmp(N);
  
  for (Subscript i = 0; i < N; i++) tmp[i] = A[i] - B[i];
  
  return tmp;
}

template <class T>
inline Vec<T> operator* (const T &val, const Vec<T> &A) 
{
  Subscript N = A.dim();
  
  Vec<T> tmp( N );
  
  for (Subscript i = 0; i < N; i++) tmp[i] = val * A[i];
  
  return tmp;
}


template<class T>
inline Vec<T> operator* (const Vec<T> &A, const T &val) 
{
  Subscript N = A.dim();
  
  Vec<T> tmp( N );
  
  for (Subscript i = 0; i < N; i++) tmp[i] = val * A[i];
  
  return tmp;
}


template<class T>
inline Vec<T> operator/ (const Vec<T> &A, const T &val) 
{
  Subscript N = A.dim();
  
  Vec<T> tmp( N );
  
  for (Subscript i = 0; i < N; i++) tmp[i] = A[i] / val;
  
  return tmp;
}

template <class T>
inline Vec<T> operator* (const Vec<T> &A, const Vec<T> &B) 
{
  Subscript N = A.dim();
  
  assert( N == B.dim() );
  
  Vec<T> tmp(N);
  
  for (Subscript i = 0; i < N; i++) tmp[i] = A[i] * B[i];
  
  return tmp;
}

template <class T>
inline T dot_prod( const Vec<T> &A, const Vec<T> &B )
{
  Subscript N = A.dim();
  assert( N == B.dim() );
  
  T sum = (T)0;
  
  for (Subscript i = 0; i < N; i++) sum += A[i] * B[i];
  
  return sum;
}


// 
// basic matrix algorithms
//

template <class T> 
inline Matrix<T> operator+ ( const Matrix<T> & A, const Matrix<T> &B )
{
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  assert( X == B.num_rows() );
  assert( Y == B.num_cols() );
  
  Matrix<T> tmp(X, Y);
  Subscript i, j;
  
  for (i = 0; i < X; i++) 
    for (j = 0; j < Y; j++)
      tmp[i][j] = A[i][j] + B[i][j];
  
  return tmp;
}

template <class T>
inline Matrix<T> operator- (const Matrix<T> &A, const Matrix<T> &B) 
{
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  assert( X == B.num_rows() );
  assert( Y == B.num_cols() );
  
  Matrix<T> tmp( X, Y );
  Subscript i, j;
  
  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++)
      tmp[i][j] = A[i][j] - B[i][j];
  
  return tmp;
}

template <class T>
inline Matrix<T> operator* (const Matrix<T> &A, const T &val) 
{
  
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  Matrix<T> tmp( X, Y );
  Subscript i, j;
  
  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++)
      tmp[i][j] = A[i][j] * val;
    
  return tmp;
}

template <class T>
inline Matrix<T> operator* (const T &val, const Matrix<T> &A) 
{
  
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  Matrix<T> tmp( X, Y );
  Subscript i, j;
  
  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++)
      tmp[i][j] = A[i][j] * val;
    
  return tmp;
}

template <class T>
inline Matrix<T> operator/ (const Matrix<T> &A, const T &val) 
{
  
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  Matrix<T> tmp( X, Y );
  Subscript i, j;
  
  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++)
      tmp[i][j] = A[i][j] / val;
    
  return tmp;
}

template <class T>
inline bool isEqual( const Matrix<T> &A, const Matrix<T> &B ) 
{
  // return true if two matrices are equal in all elements

  bool result = true;

  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  if ( (X != B.num_rows()) || (Y != B.num_cols()) ) {
    result = false;
    return result;
  }
  
  Subscript i, j;
  
  for (i = 0; i < X; i++) {
    for (j = 0; j < Y; j++) {
      if ( A[i][j] != B[i][j] ) {
	result = false;
	return result;
      }
    }
  }
    
  return result;
}

template <class T>
inline Matrix<T> mult_element( const Matrix<T> &A, const Matrix<T> &B ) 
{
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  assert( X == B.num_rows() );
  assert( Y == B.num_cols() );
  
  Matrix<T> tmp( X, Y );
  Subscript i, j;
  
  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++) 
      tmp[i][j] = A[i][j] * B[i][j];
  
  return tmp;
}

template <class T>
inline Matrix<T> div_element( const Matrix<T> &A, const Matrix<T> &B ) 
{
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  assert( X == B.num_rows() );
  assert( Y == B.num_cols() );
  
  Matrix<T> tmp( X, Y );
  Subscript i, j;
  
  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++) 
      tmp[i][j] = A[i][j] / B[i][j];
  
  return tmp;
}

template <class T>
inline void mult_element( const Matrix<T> &A, const Matrix<T> &B, Matrix<T> out ) 
{
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  assert( X == B.num_rows() );
  assert( Y == B.num_cols() );
  
  if ( (X != out.num_rows()) || (Y != out.num_cols()) ) out.newsize(X,Y);
  
  Subscript i, j;
  
  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++) 
      out[i][j] = A[i][j] * B[i][j];
}


template <class T>
inline Matrix<T> transpose( const Matrix<T> &A )
{
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  Matrix<T> tmp( Y, X );
  Subscript i, j;
  
  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++) 
      tmp[j][i] = A[i][j];
  
  return tmp;
}

template <class T>
inline void transpose( const Matrix<T> &A, Matrix<T> &out )
{
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  if ( (X != out.num_cols()) || (Y != out.num_rows()) ) out.newsize( Y, X );
    
  Subscript i, j;
  
  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++) 
      out[j][i] = A[i][j];
}

template <class T, class S>
inline void convertMatrixType( const Matrix<T> &A, Matrix<S> &out )
{
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  if ( (X != out.num_rows()) || (Y != out.num_cols()) ) out.newsize( X, Y );
    
  Subscript i, j;

  for (i = 0; i < X; i++)
    for (j = 0; j < Y; j++)
      out[i][j] = (S)A[i][j];
  
}

template <class T>
inline Matrix<T> matmult( const Matrix<T> &A, const Matrix<T> &B ) 
{
#ifdef BOUNDS_CHECK
  assert( A.num_cols() == B.num_rows() );
#endif
  
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  Subscript K = B.num_cols();
  Subscript i, j, k;
  
  Matrix<T> tmp( X, K );
  
  T sum;
  
  for (i = 0; i < X; i++) 
    for (k = 0; k < K; k++) {
      
      sum = (T)0;
      for (j = 0; j < Y; j++)
	sum += A[i][j] * B[j][k];
      
      tmp[i][k] = sum;
    }
  
  return tmp;
}

template <class T>
inline Matrix<T> operator* ( const Matrix<T> &A, const Matrix<T> &B )
{
  return matmult( A, B );
}


template <class T>
inline int matmult( Matrix<T> &C, const Matrix<T> &A, const Matrix<T> &B ) 
{
  assert( A.num_cols() == B.num_rows() );
  
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  Subscript K = B.num_cols();
  Subscript i, j, k;
  
  if ( (X != C.num_rows()) || (K != C.num_cols()) ) C.newsize( X, K );
  
  T sum;
  
  const T *row_i;
  const T *col_k;
  
  for (i = 0; i < X; i++)
    for (k = 0; k < K; k++) {
      
      row_i = &(A[i][0]);
      col_k = &(B[0][k]);
      sum = 0;
      
      for (j = 0; j < Y; j++) {
	sum += *row_i * *col_k;
	row_i++;
	col_k += K;
      }
      C[i][k] = sum;
    }
  
  return 0;
}

template <class T>
inline Vec<T> matmult( const Matrix<T> &A, const Vec<T> &x )
{
#ifdef BOUNDS_CHECK
  assert( A.num_cols() == x.dim() );
#endif
  
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  Vec<T> tmp(X);
  
  T sum;
  
  for (Subscript i = 0; i < X; i++) {
    sum = (T)0;
    
    const T *row_i = A[i];
    for (Subscript j = 0; j < Y; j++) sum += row_i[j] * x[j];
    
    tmp[i] = sum;
  }
  
  return tmp;
}

template <class T>
inline Vec<T> operator* ( const Matrix<T> &A, const Vec<T> &x )
{
  return matmult( A, x );
}

template <class T>
inline Vec<T> matmult( const Vec<T> &x, const Matrix<T> &A )
{
#ifdef BOUNDS_CHECK
  assert( A.num_rows() == x.dim() );
#endif
  
  Subscript X = A.num_rows();
  Subscript Y = A.num_cols();
  
  Vec<T> tmp(X);
  
  T sum;
  
  for (Subscript i = 0; i < Y; i++) {
    sum = (T)0;
    
    for (Subscript j = 0; j < X; j++) sum += A[j][i] * x[j];
    
    tmp[i] = sum;
  }
  
  return tmp;
}

template <class T>
inline Vec<T> operator* ( const Vec<T> &x, const Matrix<T> &A )
{
  return matmult( x, A );
}

template <class T>
inline Matrix<T> eye( Subscript n, T a )
{
  // create a diagonal matrix of size nxn who's diagonal
  // entries all have value "a"
  Matrix<T> p( n, n, (T)0 );
  
  for (Subscript i = 0; i < n; i++) p[i][i] = a;
  
  return p;
}

template <class T>
inline Matrix<T> eye( const Vec<T> &x ) 
{
  // create a diagonal matrix who's diagonal is equal to vector x
  Subscript n = x.dim();
  
  Matrix<T> p( n, n, (T)0 );
  
  for (Subscript i = 0; i < n; i++) p[i][i] = x[i];
  
  return p;
}

template <class T>
inline Matrix<T> eye( const Vec<T> &x, int k ) 
{
  // create a diagonal matrix who's diagonal is equal to vector x
  Subscript n = x.dim();
  int off = abs(k);
  
  Matrix<T> p( n+off, n+off, (T)0 );
  
  if ( k < 0 ) {
    for (Subscript i = 0; i < n; i++) p[i+off][i] = x[i];
  }
  else {
    for (Subscript i = 0; i < n; i++) p[i][i+off] = x[i];
  }
  
  return p;
}

template <class T>
inline Matrix<T> eye( const Vec<T> &x, int k, int beginD, int endD ) 
{
  // create a diagonal matrix who's diagonal is equal to vector x
  Subscript n = (Subscript) (endD - beginD + 1);
  int off = k;
  if ( off < 0 ) off = -k;
  
  Matrix<T> p( n+off, n+off, (T)0 );
  
  if ( k < 0 ) {
    for (Subscript i = 0; i < n; i++) p[i+off][i] = x[beginD + i];
  }
  else {
    for (Subscript i = 0; i < n; i++) p[i][i+off] = x[beginD + i];
  }
  
  return p;
}

template < class T >
inline Vec<T> lookUp( const Matrix<T> &z,
		      const Vec<T> &xi, const Vec<T> &yi ) 
{
  //
  // this is an replacement for matlab's interp2
  // linear interpolation
  //
  Subscript xi_x = xi.dim();
  Subscript yi_x = yi.dim();
  Subscript z_x = z.num_rows();
  Subscript z_y = z.num_cols();
  Subscript indx, indy;
  T u, v;
  
  assert( xi_x == yi_x );
  
  Vec<T> tmp( xi_x );
  
  for (Subscript i = 0; i < xi_x; i++) {

#ifdef BOUNDS_CHECK
    assert( ( xi[i] >= (T)0 ) && ( xi[i] <= (T)(z_y-1) ) );
    assert( ( yi[i] >= (T)0 ) && ( yi[i] <= (T)(z_x-1) ) );
#endif
  
    indx = (int) floor(yi[i]);
    indy = (int) floor(xi[i]);
  
    u = z[ indx ][ indy ] +
      ( z[ indx + 1 ][ indy ] - z[ indx ][ indy ] ) * ( yi[i] - (T)indx );
    v = z[ indx ][ indy + 1 ] +
      ( z[ indx + 1 ][ indy + 1 ] - z[ indx ][ indy + 1 ] ) *
      ( yi[i] - (T)indx );
      
    tmp[i] = u + ( v - u ) * ( xi[i] - (T)indy );
  }
  
  return tmp;
}

template < class T >
inline Vec<T> lookUp( const Matrix<T> &z,
		      const Vec<Subscript> &xi, const Vec<Subscript> &yi ) 
{
  //
  // this is an replacement for matlab's interp2
  // simply a lookup table, without interpolation
  //
  Subscript xi_x = xi.dim();
  Subscript yi_x = yi.dim();
  Subscript z_x = z.num_rows();
  Subscript z_y = z.num_cols();
  
  assert( xi_x == yi_x );
  
  Vec<T> tmp( xi_x );
  
  for (Subscript i = 0; i < xi_x; i++) {

#ifdef BOUNDS_CHECK
    assert( ( xi[i] >= 0 ) && ( xi[i] < z_x ) );
    assert( ( yi[i] >= 0 ) && ( yi[i] < z_y ) );
#endif
  
    tmp[i] = z[ xi[i] ][ yi[i] ];
  }
  
  return tmp;
}


inline Matrix<double> rotMatYZX( double gamma, double beta, double alpha )
{
  // gamma, rotation about Y axis (flexion)
  // beta, rotation about Z axis (internal-external angulation)
  // alpha, rotation about X axis (varus-valgus angulation)
  //
  // all inputs are in degree
  double pii = atan(1.0) * 4.0;
  double a = alpha * pii / 180.0, ca = cos(a), sa = sin(a);
  double b =  beta * pii / 180.0, cb = cos(b), sb = sin(b);
  double g = gamma * pii / 180.0, cg = cos(g), sg = sin(g);
  Matrix<double> tmp( 3, 3 );

  // first row
  tmp[0][0] = cb * cg;
  tmp[0][1] = -sb;
  tmp[0][2] = cb * sg;
  
  // second row
  tmp[1][0] = ca * sb * cg + sa * sg;
  tmp[1][1] = ca * cb;
  tmp[1][2] = ca * sb * sg - sa * cg;
  
  // third row
  tmp[2][0] = sa * sb * cg - ca * sg;
  tmp[2][1] = sa * cb;
  tmp[2][2] = sa * sb * sg + ca * cg;
    
  return (tmp);
}

inline Matrix<double> rotMatX( double theta )
{
  // all inputs are in degree
  double pii = atan(1.0) * 4.0;
  double t = theta * pii / 180.0, ct = cos(t), st = sin(t);
  
  Matrix<double> tmp( 3, 3 );

  tmp[0][0] = 1.0;  tmp[0][1] = 0.0; tmp[0][2] = 0.0;
  tmp[1][0] = 0.0;  tmp[1][1] = ct;  tmp[1][2] = -st;
  tmp[2][0] = 0.0;  tmp[2][1] = st;  tmp[2][2] = ct;
      
  return (tmp);
}

inline Matrix<double> rotMatY( double theta )
{
  // all inputs are in degree
  double pii = atan(1.0) * 4.0;
  double t = theta * pii / 180.0, ct = cos(t), st = sin(t);
  
  Matrix<double> tmp( 3, 3 );

  tmp[0][0] = ct;   tmp[0][1] = 0.0;  tmp[0][2] = st;
  tmp[1][0] = 0.0;  tmp[1][1] = 1.0;  tmp[1][2] = 0.0;
  tmp[2][0] = -st;  tmp[2][1] = 0.0;  tmp[2][2] = ct;
      
  return (tmp);
}

inline Matrix<double> rotMatZ( double theta )
{
  // all inputs are in degree
  double pii = atan(1.0) * 4.0;
  double t = theta * pii / 180.0, ct = cos(t), st = sin(t);
  
  Matrix<double> tmp( 3, 3 );

  tmp[0][0] = ct;   tmp[0][1] = -st;  tmp[0][2] = 0.0;
  tmp[1][0] = st;   tmp[1][1] = ct;   tmp[1][2] = 0.0;
  tmp[2][0] = 0.0;  tmp[2][1] = 0.0;  tmp[2][2] = 1.0;
      
  return (tmp);
}

#endif
