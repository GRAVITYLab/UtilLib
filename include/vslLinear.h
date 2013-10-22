/** 
 *  vslLinear.h
 *  This file defines interface to graphics linear operation classes for VSL
 *  
 *  Author: Daqing Xue
 *  Date: 7/16/2005
 */

/**
 * Interface to graphics linear operation classes
 * 
 *   Daqing Xue
 *   Department of Computer & Information Science
 *   The Ohio State University
 *   Apr. 2003
 *
 * Copyright (C) 2003 - 2005 Daqing Xue
 *
 * The source code can be freely used, redistributed, and modified for 
 * any purpose, provided that the above copyright information appears 
 * in all copies.
 */

#ifndef __VSL_LINEAR_H
#define __VSL_LINEAR_H

#include <cmath>
#include <iostream>

namespace utillib
{

#pragma warning( push )
#pragma warning( disable : 4201 ) // Nonstandard extension used : nameless struct/union

#define VSL_EPSILON_AS_ZERO 0.000005

template <typename T> inline
  bool is_zero ( T t, T epsilon = VSL_EPSILON_AS_ZERO )
{
  return fabs(t) < epsilon;
}


/**
 * 
 *  2D vector template class
 *
 */

template <typename T>
class vec2
{
 public:
  union {
    T _v[2];
    struct {
      T x, y;
    };
    struct {
      T u, v;
    };
    struct {
      T s, t;
    };
  };

 public:
  explicit vec2( const T & t = T() ) 
    {  _v[0] = t; _v[1] = t;  }
  
  explicit vec2( const T * t )
    {  _v[0] = t[0]; _v[1] = t[1];  }
  
  vec2( const vec2<T> & v ) 
    {  _v[0] = v[0]; _v[1] = v[1];  }

  vec2( const T & x, const T & y )  
    {  _v[0] = x; _v[1] = y;  }

  const T * getValue( T * t = NULL ) const
    {  if (t!=NULL) { t[0] = _v[0]; t[1] = _v[1];  }   return _v; }
  
  void setValue( const T * t )
    {  _v[0] = t[0]; _v[1] = t[1];  }

  void getValue( T & x, T & y ) const
    {  x = _v[0]; y = _v[1];  }
  
  void setValue( const T & x, const T & y )
    {  _v[0] = x; _v[1] = y;  }

  void reset( const T & t )
    {  _v[0] = t; _v[1] = t;  }

  void zero()
    {  _v[0] = T(0); _v[1] = T(0);  }
  
  void negate()
    {  _v[0] = -_v[0]; _v[1] = -_v[1];  }
  
  bool isZero() const
    {  return is_zero(_v[0]) && is_zero(_v[1]);  }

  bool agree( const vec2<T> & v) const
    {  return (*this - v).isZero();  }

  T length() const
    {  return T(sqrt(length2()));  }	
  
  T length2() const
    {  return T(_v[0]*_v[0] + _v[1]*_v[1]);  }	

  T distance(const vec2<T> & v) const
    { return (*this-v).length();  }

  T dot( const vec2<T> & v ) const
    {  return T(_v[0]*v[0] + _v[1]*v[1]);  }
  
  void normalize() 
    { 
      T l;
      l = length();
      if (l != T(0)) { T rl = T(1)/l; _v[0] *= rl; _v[1] *= rl;  }
    }
  
  void print() const
    {  std::cout << "[ " << _v[0] << " " << _v[1] << " ]\n";  }
  
  // operator functions
  
  T operator () ( int i ) const
    { return _v[i]; }
  
  T & operator [] ( int i )
    { return _v[i]; }
  
  const T & operator [] ( int i ) const
    { return _v[i]; }

  vec2<T> & operator = ( const vec2<T> & v )
    { _v[0] = v[0];  _v[1] = v[1];    return (*this); }
  
  vec2<T> & operator += ( const vec2<T> & v )
    { _v[0] += v[0];  _v[1] += v[1];  return (*this); }
  
  vec2<T> & operator -= ( const vec2<T> & v )
    { _v[0] -= v[0];  _v[1] -= v[1];  return (*this); }
  
  vec2<T> & operator *= ( const vec2<T> & v )
    { _v[0] *= v[0];  _v[1] *= v[1];  return (*this); }

  vec2<T> & operator *= ( T t )
    { _v[0] *= t;  _v[1] *= t;  return (*this); }

  vec2<T> & operator /= ( T t )
    { if (t != T(0)) { T rt = T(1)/t;  _v[0] *= rt; _v[1] *= rt; }  return (*this);  }
  
  vec2<T> operator + () const
    { return vec2<T>(*this); }
  
  vec2<T> operator - () const
    { vec2<T> rv(*this); rv.negate(); return rv; }
  
  vec2<T> operator + ( const vec2<T> &v ) const
    { vec2<T> rv(*this); return rv += v; }
  
  vec2<T> operator - ( const vec2<T> &v ) const
    { vec2<T> rv(*this); return rv -= v; }

  vec2<T> operator * ( const vec2<T> &v ) const
    { vec2<T> rv(*this); return rv *= v; }

  vec2<T> operator * ( T t ) const
    { vec2<T> rv(*this); return rv *= t; }
  
  vec2<T> operator / ( T t ) const
    { vec2<T> rv(*this); return rv /= t; }

  bool operator == ( const vec2<T> & v ) const
    {  return _v[0]==v[0] && _v[1]==v[1]; }

  bool operator != ( const vec2<T> & v) const
    {  return !(*this == v);  }
};

// vec2<T> operator functions
		
template <typename T> inline
  vec2<T> operator * ( T t, const vec2<T> & v )
{ 
  return v*t; 
}

template <typename T> inline
  T dot ( const vec2<T> & v, const vec2<T> & u )
{
  return v.dot(u);
}

template <typename T> inline
  vec2<T> normalize ( const vec2<T> & v )
{
  vec2<T> rv(v);
  rv.normalize();
  return rv;
}

/**
 * 
 *  3D vector template class
 *
 */

template <typename T>
class vec3
{
 public:
  union {
    T _v[3];
    struct {
      T x, y, z;
    };
    struct {
      T r, g, b;
    };
  };

 public:
  explicit vec3( const T & t = T() ) 
    {  _v[0] = t; _v[1] = t;  _v[2] = t;  }
  
  explicit vec3( const T * t )
    {  _v[0] = t[0]; _v[1] = t[1];  _v[2] = t[2];  }

  explicit vec3( const vec2<T> & v, const T t = T() ) 
    {  _v[0] = v[0]; _v[1] = v[1]; _v[2] = t;  }

  vec3( const T & x, const T & y, const T & z )  
    {  _v[0] = x; _v[1] = y;  _v[2] = z;  }

  vec3( const vec3<T> & v ) 
    {  _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2];  }

  const T * getValue( T * t = NULL ) const
    {  if (t!=NULL) { t[0] = _v[0]; t[1] = _v[1]; t[2] = _v[2]; }   return _v; }
  
  void setValue( const T * t )
    {  _v[0] = t[0]; _v[1] = t[1]; _v[2] = t[2];  }

  void getValue( T & x, T & y, T & z ) const
    {  x = _v[0]; y = _v[1]; z = _v[2];  }
  
  void setValue( const T & x, const T & y, const T & z )
    {  _v[0] = x; _v[1] = y; _v[2] = z;  }

  void reset( const T & t )
    {  _v[0] = t; _v[1] = t; _v[2] = t;  }

  void zero()
    {  _v[0] = T(0); _v[1] = T(0); _v[2] = T(0); }
  
  void negate()
    {  _v[0] = -_v[0]; _v[1] = -_v[1]; _v[2] = -_v[2];  }
  
  bool isZero() const
    {  return is_zero(_v[0]) && is_zero(_v[1]) && is_zero(_v[2]);  }

  bool agree( const vec3<T> & v) const
    {  return (*this - v).isZero();  }

  T length() const
    {  return T(sqrt(length2()));  }	
  
  T length2() const
    {  return T(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2]);  }	

  T distance(const vec3<T> & v) const
    { return (*this-v).length();  }

  T dot( const vec3<T> & v ) const
    {  return T(_v[0]*v[0] + _v[1]*v[1] + _v[2]*v[2]);  }
  
  void normalize() 
    { 
      T l;
      l = length();
      if (l != T(0)) { T rl = T(1)/l; _v[0] *= rl; _v[1] *= rl; _v[2] *= rl; }
    }

  vec3<T> cross(const vec3<T> & rhs ) const
    {
      vec3<T> rv;
      rv[0] = _v[1]*rhs[2] - _v[2]*rhs[1];
      rv[1] = _v[2]*rhs[0] - _v[0]*rhs[2];
      rv[2] = _v[0]*rhs[1] - _v[1]*rhs[0];	
      return rv;
    }

  void print() const
    {  std::cout << "[ " << _v[0] << " " << _v[1] <<" " << _v[2] << " ]\n";  }
  
  // operator functions
  
  T operator () ( int i ) const
    { return _v[i]; }
  
  T & operator [] ( int i )
    { return _v[i]; }
  
  const T & operator [] ( int i ) const
    { return _v[i]; }

  vec3<T> & operator = ( const vec3<T> & v )
    { _v[0] = v[0];  _v[1] = v[1];  _v[2] = v[2];     return (*this); }
  
  vec3<T> & operator += ( const vec3<T> & v )
    { _v[0] += v[0];  _v[1] += v[1];  _v[2] += v[2];  return (*this); }
  
  vec3<T> & operator -= ( const vec3<T> & v )
    { _v[0] -= v[0];  _v[1] -= v[1];  _v[2] -= v[2];  return (*this); }
  
  vec3<T> & operator *= ( const vec3<T> & v )
    { _v[0] *= v[0];  _v[1] *= v[1];  _v[2] *= v[2];  return (*this); }

  vec3<T> & operator *= ( T t )
    { _v[0] *= t;  _v[1] *= t;  _v[2] *= t;  return (*this); }
  
  vec3<T> & operator /= ( T t )
    {
      if (t != T(0)) { T rt = T(1)/t; _v[0] *= rt; _v[1] *= rt; _v[2] *= rt; }
      return (*this);
    }
  
  vec3<T> operator + () const
    { return vec3<T>(*this); }
  
  vec3<T> operator - () const
    { vec3<T> rv(*this); rv.negate(); return rv; }
  
  vec3<T> operator + ( const vec3<T> &v ) const
    { vec3<T> rv(*this); return rv += v; }
  
  vec3<T> operator - ( const vec3<T> &v ) const
    { vec3<T> rv(*this); return rv -= v; }

  vec3<T> operator * ( const vec3<T> &v ) const
    { vec3<T> rv(*this); return rv *= v; }

  vec3<T> operator * ( T t ) const
    { vec3<T> rv(*this); return rv *= t; }
  
  vec3<T> operator / ( T t ) const
    { vec3<T> rv(*this); return rv /= t; }

  bool operator == ( const vec3<T> & v ) const
    {  return _v[0]==v[0] && _v[1]==v[1] && _v[2]==v[2]; }

  bool operator != ( const vec3<T> & v) const
    {  return !(*this == v);  }
};

// vec3<T> operator functions
		
template <typename T> inline
  vec3<T> operator * ( T t, const vec3<T> & v )
{ 
  return v*t; 
}
template <typename T> inline
  T dot ( const vec3<T> & v, const vec3<T> & u )
{
  return v.dot(u);
}

template <typename T> inline
  vec3<T> normalize ( const vec3<T> & v )
{
  vec3<T> rv(v);
  rv.normalize();
  return rv;
}

template <typename T> inline
  vec3<T> cross ( const vec3<T> & lhs, const vec3<T> rhs )
{
  return lhs.cross(rhs);
}

/**
 * 
 *  4D vector template class
 *
 */

template <typename T>
class vec4
{
 public:
  union {
    T _v[4];
    struct {
      T x, y, z, w;
    };
    struct {
      T r, g, b, a;
    };
  };

 public:
  explicit vec4( const T & t = T() ) 
    {  _v[0] = t; _v[1] = t; _v[2] = t; _v[3] = t;  }
  
  explicit vec4( const T * t )
    {  _v[0] = t[0]; _v[1] = t[1]; _v[2] = t[2]; _v[3] = t[3];  }

  explicit vec4( const vec3<T> & v, const T t = T() ) 
    {  _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2]; _v[3] = t;  }

  vec4( const T & x, const T & y, const T & z, const T & w )  
    {  _v[0] = x; _v[1] = y; _v[2] = z; _v[3] = w;  }

  vec4( const vec4<T> & v ) 
    {  _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2]; _v[3] = v[3];  }

  const T * getValue( T * t = NULL ) const
    {  
      if (t!=NULL) { t[0] = _v[0]; t[1] = _v[1]; t[2] = _v[2]; t[3] = _v[3]; }
      return _v; 
    }
  
  void setValue( const T * t )
    {  _v[0] = t[0]; _v[1] = t[1]; _v[2] = t[2]; _v[3] = t[3];  }

  void getValue( T & x, T & y, T & z, T & w ) const
    {  x = _v[0]; y = _v[1]; z = _v[2]; w = _v[3];  }
  
  void setValue( const T & x, const T & y, const T & z, const T & w )
    {  _v[0] = x; _v[1] = y; _v[2] = z; _v[3] = w;  }

  void reset( const T & t )
    {  _v[0] = t; _v[1] = t; _v[2] = t; _v[3] = t;  }

  void zero()
    {  _v[0] = T(0); _v[1] = T(0); _v[2] = T(0); _v[3] = T(0);  }
  
  void negate()
    {  _v[0] = -_v[0]; _v[1] = -_v[1]; _v[2] = -_v[2]; _v[3] = -_v[3];  }
  
  bool isZero() const
    {  return is_zero(_v[0]) && is_zero(_v[1]) && is_zero(_v[2]) && is_zero(_v[3]);  }

  bool agree( const vec4<T> & v) const
    {  return (*this - v).isZero();  }
  
  T length() const
    {  return T(sqrt(length2()));  }	
  
  T length2() const
    {  return T(_v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]);  }	

  T distance(const vec4<T> & v) const
    { return (*this-v).length();  }

  T dot( const vec4<T> & v ) const
    {  return T(_v[0]*v[0] + _v[1]*v[1] + _v[2]*v[2] + _v[3]*v[3]);  }
  
  void normalize() 
    { 
      T l;
      l = length();
      if (l != T(0)) { T rl = T(1)/l; _v[0] *= rl; _v[1] *= rl; _v[2] *= rl; _v[3] *= rl; }
    }
  void homogenize()
    {  
      if (_v[3] != T(0)) { T r = T(1)/_v[3]; _v[0] *= r; _v[1] *= r; _v[2] *= r; }
      _v[3] = 1.0; 
    }

  void print() const
    {  std::cout << "[ " << _v[0] << " " << _v[1] <<" " << _v[2] << " " << _v[3] << " ]\n";  }
  
  // operator functions
  
  T operator () ( int i ) const
    { return _v[i]; }
  
  T & operator [] ( int i )
    { return _v[i]; }
  
  const T & operator [] ( int i ) const
    { return _v[i]; }

  vec4<T> & operator = ( const vec4<T> & v )
    { _v[0] = v[0];  _v[1] = v[1];  _v[2] = v[2];  _v[3] = v[3];      return (*this); }
  
  vec4<T> & operator += ( const vec4<T> & v )
    { _v[0] += v[0];  _v[1] += v[1];  _v[2] += v[2];  _v[3] += v[3];  return (*this); }
  
  vec4<T> & operator -= ( const vec4<T> & v )
    { _v[0] -= v[0];  _v[1] -= v[1];  _v[2] -= v[2];  _v[3] -= v[3];  return (*this); }
  
  vec4<T> & operator *= ( const vec4<T> & v )
    { _v[0] *= v[0];  _v[1] *= v[1];  _v[2] *= v[2];  _v[3] *= v[3];  return (*this); }

  vec4<T> & operator *= ( T t )
    { _v[0] *= t;  _v[1] *= t;  _v[2] *= t;  _v[3] *= t;  return (*this); }
  
  vec4<T> & operator /= ( T t )
    {
      if (t != T(0)) { T rt = T(1)/t; _v[0] *= rt; _v[1] *= rt; _v[2] *= rt; _v[3] *= rt; }
      return (*this);
    }
  
  vec4<T> operator + () const
    { return vec4<T>(*this); }
  
  vec4<T> operator - () const
    { vec4<T> rv(*this); rv.negate(); return rv; }
  
  vec4<T> operator + ( const vec4<T> &v ) const
    { vec4<T> rv(*this); return rv += v; }
  
  vec4<T> operator - ( const vec4<T> &v ) const
    { vec4<T> rv(*this); return rv -= v; }

  vec4<T> operator * ( const vec4<T> &v ) const
    { vec4<T> rv(*this); return rv *= v; }

  vec4<T> operator * ( T t ) const
    { vec4<T> rv(*this); return rv *= t; }
  
  vec4<T> operator / ( T t ) const
    { vec4<T> rv(*this); return rv /= t; }

  bool operator == ( const vec4<T> & v ) const
    {  return _v[0]==v[0] && _v[1]==v[1] && _v[2]==v[2] && _v[3]==v[3]; }

  bool operator != ( const vec4<T> & v) const
    {  return !(*this == v);  }
};

// vec4<T> operator functions
		
template <typename T> inline
  vec4<T> operator * ( T t, const vec4<T> & v )
{ 
  return v*t; 
}
template <typename T> inline
  T dot ( const vec4<T> & v, const vec4<T> & u )
{
  return v.dot(u);
}

template <typename T> inline
  vec4<T> normalize ( const vec4<T> & v )
{
  vec4<T> rv(v);
  rv.normalize();
  return rv;
}

template <typename T> inline
  vec4<T> cross ( const vec4<T> & lhs, const vec4<T> rhs )
{
  return lhs.cross(rhs);
}

template <typename T> inline
  vec3<T> homogenize ( const vec4<T> & v )
{
  vec4<T> rv(v);
  rv.homogenize();
  return vec3<T>(&rv[0]);
}

/**
 * 
 *  3D matrix template class
 *
 */

template <typename T>
class mat3
{
 protected:
  union {
    // 1D array observes OpenGL matrix element order, aka, column-major.
    T _m[9];
    struct {
      T a00, a10, a20,             // [ a00(0) a01(3) a02(6) ] //
        a01, a11, a21,             // [ a10(1) a11(4) a12(7) ] //
        a02, a12, a22;             // [ a20(2) a21(5) a22(8) ] //
    };
  };

 public:
  explicit mat3() 
    { identity();  }

  explicit mat3( const T & t ) 
    { reset( t );     }
  
  explicit mat3( const T * t )
    { setValue ( t ); }

  mat3( const T & m00, const T & m01, const T & m02,   
        const T & m10, const T & m11, const T & m12,   
        const T & m20, const T & m21, const T & m22 )  
    {
      setValue( m00, m01, m02,
                m10, m11, m12,  
                m20, m21, m22 );
    }

  mat3( const mat3<T> & m ) 
    {
      _m[0] = m[0]; _m[3] = m[3];  _m[6] = m[6];  
      _m[1] = m[1]; _m[4] = m[4];  _m[7] = m[7];  
      _m[2] = m[2]; _m[5] = m[5];  _m[8] = m[8];  
    }

  mat3( const vec3<T> & v0, const vec3<T> & v1, const vec3<T> & v2  ) 
    {
      _m[0] = v0[0]; _m[3] = v0[1];  _m[6] = v0[2];  
      _m[1] = v1[0]; _m[4] = v1[1];  _m[7] = v1[2];  
      _m[2] = v2[0]; _m[5] = v2[1];  _m[8] = v2[2];  
    }

  const T * getValue( T * t = NULL ) const
    // return 1D array observing OpenGL matrix element order, aka, column-major
    {  
      if (t!=NULL) { 
        t[0] = _m[0]; t[3] = _m[3]; t[6] = _m[6]; 
        t[1] = _m[1]; t[4] = _m[4]; t[7] = _m[7]; 
        t[2] = _m[2]; t[5] = _m[5]; t[8] = _m[8]; 
      }   
      return _m; 
    }

  void setValue( const T * t )
    {  
      _m[0] = t[0]; _m[3] = t[3]; _m[6] = t[6];  
      _m[1] = t[1]; _m[4] = t[4]; _m[7] = t[7];  
      _m[2] = t[2]; _m[5] = t[5]; _m[8] = t[8];  
    }
  
  void getValue(  T & m11, T & m12, T & m13,   
                  T & m21, T & m22, T & m23,   
                  T & m31, T & m32, T & m33 ) const
    {
      m11 = _m[0]; m12 = _m[3];  m13 = _m[6];  
      m21 = _m[1]; m22 = _m[4];  m23 = _m[7];  
      m31 = _m[2]; m32 = _m[5];  m33 = _m[8];  
    }

  void setValue(  const T & m00, const T & m01, const T & m02,   
                  const T & m10, const T & m11, const T & m12,   
                  const T & m20, const T & m21, const T & m22 )  
    {
      _m[0] = m00; _m[3] = m01;  _m[6] = m02;  
      _m[1] = m10; _m[4] = m11;  _m[7] = m12;  
      _m[2] = m20; _m[5] = m21;  _m[8] = m22;  
    }

  vec3<T> getRowVector( int i ) const
    { return vec3<T>(_m[i], _m[i+3], _m[i+6]); }

  void setRowVector( int i, const vec3<T> & v )
    { _m[i] = v[0]; _m[i+3] = v[1]; _m[i+6] = v[2]; }

  void getRowValue( int i , T * t ) const
    { t[0] = _m[i]; t[1] = _m[i+3]; t[2] = _m[i+6]; }

  void setRowValue( int i , const T * t )
    { _m[i] = t[0]; _m[i+3] = t[1]; _m[i+6] = t[2]; }
    
  vec3<T> getColVector( int j ) const
    { return vec3<T>(&_m[j*3]); }

  void setColVector( int j, const vec3<T> & v )
    { _m[j*3] = v[0]; _m[j*3+1] = v[1]; _m[j*3+2] = v[2]; }

  void getColValue( int j , T * t ) const
    { t[0] = _m[j*3]; t[1] = _m[j*3+1]; t[2] = _m[j*3+2]; }

  void setColValue( int j, const T * t )
    { _m[j*3] = t[0]; _m[j*3+1] = t[1]; _m[j*3+2] = t[2]; }

  void reset( const T & t )
  {
      _m[0] = t; _m[3] = t;  _m[6] = t;  
      _m[1] = t; _m[4] = t;  _m[7] = t;  
      _m[2] = t; _m[5] = t;  _m[8] = t;  
  }

  void zero()
    {  reset( 0 ); }
  
  void negate()
    {
      _m[0] = -_m[0]; _m[3] = -_m[3];  _m[6] = -_m[6];  
      _m[1] = -_m[1]; _m[4] = -_m[4];  _m[7] = -_m[7];  
      _m[2] = -_m[2]; _m[5] = -_m[5];  _m[8] = -_m[8];  
    }
  
  bool isZero() const
    {  
      return  is_zero(_m[0]) && is_zero(_m[1]) && is_zero(_m[2]) && 
              is_zero(_m[3]) && is_zero(_m[4]) && is_zero(_m[5]) && 
              is_zero(_m[6]) && is_zero(_m[7]) && is_zero(_m[8]);
    }

  bool agree( const mat3<T> & m) const
    {  return (*this - m).isZero();  }

  void identity()
    { 
      _m[0] = T(1); _m[3] = 0;     _m[6] = 0;  
      _m[1] = 0;    _m[4] = T(1);  _m[7] = 0;  
      _m[2] = 0;    _m[5] = 0;     _m[8] = T(1);  
    }
  
  void transpose()
    { 
      T t;
      t = a01; a01 = a10; a10 = t;
      t = a12; a12 = a21; a21 = t;
      t = a20; a20 = a02; a02 = t;
    }
/*
  T minor(int r0, int r1, int c0, int c1) const
    // calculate the minor(r,c), where r = {0,1,2}-{r0,r1}, c = {0,1,2}-{c0,c1}
    { return (*this)(r0,c0) * (*this)(r1,c1) - (*this)(r0,c1) * (*this)(r1,c0); }

  T determinant() const
    { return _m[0]*minor(1,2, 1,2) - _m[1]*minor(0,2, 1,2) + _m[2]*minor(0,1, 1,2); }

  void adjoint()
    {
      mat3<T> m(*this);
#define M(r0,r1, c0,c1)  m.minor(r0,r1, c0,c1)
      _m[0] =   M(1,2, 1,2); _m[3] = - M(0,2, 1,2); _m[6] =   M(0,1, 1,2);
      _m[1] = - M(1,2, 0,2); _m[4] =   M(0,2, 0,2); _m[7] = - M(0,1, 0,2);
      _m[2] =   M(1,2, 0,1); _m[5] = - M(0,2, 0,1); _m[8] =   M(0,1, 0,1);
#undef M
    }
*/
  bool inverse()
    {
      /*
	  T det = determinant();

      if (!is_zero(det)) {
        this->adjoint();
        *this *= T(1)/det;
        return true;
      }
      */
      return false;
    }

  bool isSingular () const
    {
     /*
	  T det = determinant();
      return is_zero( det );
      */
	  return false;
    }

  T rowDotCol(int r, mat3<T> rhs, int c) const
    { 
      return (*this)(r,0) * rhs(0,c) + (*this)(r,1) * rhs(1,c) + (*this)(r,2) * rhs(2,c); 
    }

  void print() const
    {  
      std::cout << "[ " << _m[0] << " " << _m[3] <<" " << _m[6] << " ]\n"
                << "[ " << _m[1] << " " << _m[4] <<" " << _m[7] << " ]\n"
                << "[ " << _m[2] << " " << _m[5] <<" " << _m[8] << " ]\n";  
    }
  
  // operator functions
  
  T & operator () ( int r, int c )
    { return _m[c*3+r]; }

  const T & operator () ( int r, int c ) const
    { return _m[c*3+r]; }

  T & operator [] ( int i )
    { return _m[i]; }
  
  const T & operator [] ( int i ) const
    { return _m[i]; }

  mat3<T> & operator = ( const mat3<T> & m )
    { 
      _m[0] = m[0]; _m[3] = m[3];  _m[6] = m[6];  
      _m[1] = m[1]; _m[4] = m[4];  _m[7] = m[7];  
      _m[2] = m[2]; _m[5] = m[5];  _m[8] = m[8];  
      return (*this); 
    }
  
  mat3<T> & operator += ( const mat3<T> & m )
    { 
      _m[0] += m[0]; _m[3] += m[3];  _m[6] += m[6];  
      _m[1] += m[1]; _m[4] += m[4];  _m[7] += m[7];  
      _m[2] += m[2]; _m[5] += m[5];  _m[8] += m[8];  
      return (*this); 
    }
  
  mat3<T> & operator -= ( const mat3<T> & m )
    { 
      _m[0] -= m[0]; _m[3] -= m[3];  _m[6] -= m[6];  
      _m[1] -= m[1]; _m[4] -= m[4];  _m[7] -= m[7];  
      _m[2] -= m[2]; _m[5] -= m[5];  _m[8] -= m[8];  
      return (*this); 
    }

  mat3<T> & operator *= ( const mat3<T> & rhs )
    { 
      mat3<T> lhs(*this);
#define D(r, c)  lhs.rowDotCol(r, rhs, c)
      _m[0] = D(0, 0); _m[3] = D(0, 1); _m[6] = D(0, 2); 
      _m[1] = D(1, 0); _m[4] = D(1, 1); _m[7] = D(1, 2); 
      _m[2] = D(2, 0); _m[5] = D(2, 1); _m[8] = D(2, 2); 
#undef D
      return (*this); 
    }

  mat3<T> & operator *= ( T t )
    { 
      _m[0] *= t;  _m[3] *= t;  _m[6] *= t;  
      _m[1] *= t;  _m[4] *= t;  _m[7] *= t;  
      _m[2] *= t;  _m[5] *= t;  _m[8] *= t;  
      return (*this); 
    }
  
  mat3<T> & operator /= ( T t )
    { 
      if (t != T(0)) {
        T rt = T(1)/t;
        _m[0] *= rt;  _m[3] *= rt;  _m[6] *= rt;  
        _m[1] *= rt;  _m[4] *= rt;  _m[7] *= rt;  
        _m[2] *= rt;  _m[5] *= rt;  _m[8] *= rt;  
      }
      return (*this); 
    }
  
  mat3<T> operator + () const
    { return mat3<T>(*this); }
  
  mat3<T> operator - () const
    { mat3<T> rm(*this); rm.negate(); return rm; }
  
  mat3<T> operator + ( const mat3<T> & m ) const
    { mat3<T> rm(*this); return rm += m; }
  
  mat3<T> operator - ( const mat3<T> & m ) const
    { mat3<T> rm(*this); return rm -= m; }

  mat3<T> operator * ( const mat3<T> & m ) const
    { mat3<T> rm(*this); return rm *= m; }

  vec3<T> operator * ( const vec3<T> & v ) const
    { 
      return vec3<T>( getRowVector(0).dot(v), 
                      getRowVector(1).dot(v), 
                      getRowVector(2).dot(v) ); 
    }

  mat3<T> operator * ( T t ) const
    { mat3<T> rm(*this); return rm *= t; }
  
  mat3<T> operator / ( T t ) const
    { mat3<T> rm(*this); return rm /= t; }

  bool operator == ( const mat3<T> & m ) const
    {  
      return  _m[0]==m[0] && _m[3]==m[3] && _m[6]==m[6] && 
              _m[1]==m[1] && _m[4]==m[4] && _m[7]==m[7] && 
              _m[2]==m[2] && _m[5]==m[5] && _m[8]==m[8] ; 
    }

  bool operator != ( const mat3<T> & m) const
    {  return !(*this == m);  }
};

// mat3<T> operator functions
		
template <typename T> inline
  mat3<T> operator * ( T t, const mat3<T> & m )
{ 
  return m*t; 
}

template <typename T> inline
  vec3<T> operator * ( const vec3<T> & v, const mat3<T> & m )
{ 
  return vec3<T>( v.dot(m.getColVector(0)), 
                  v.dot(m.getColVector(1)), 
                  v.dot(m.getColVector(2)) ); 
}

// mat3<T> template functions

template <typename T> inline
  mat3<T>  transpose( const mat3<T> & m )  
{
  mat3<T> rm(m);
  rm.transpose();
  return rm;
}

template <typename T> inline
  mat3<T>  inverse( const mat3<T> & m )  
{
  mat3<T> rm(m);
  rm.inverse();
  return rm;
}


/**
 * 
 *  4D matrix template class
 *
 */

template <typename T>
class mat4
{
 protected:
  union {
    // 1D array observes OpenGL matrix element order, aka, column-major.
    T _m[16];
    struct {
      T a00, a10, a20, a30,             // [ a00(0) a01(4) a02(8)  a03(12) ] // 
        a01, a11, a21, a31,             // [ a10(1) a11(5) a12(9)  a13(13) ] // 
        a02, a12, a22, a32,             // [ a20(2) a21(6) a22(10) a23(14) ] // 
        a03, a13, a23, a33;             // [ a30(3) a31(7) a32(11) a33(15) ] //   
    };
  };

 public:
  explicit mat4() 
    { identity();  }

  explicit mat4( const T & t ) 
    { reset( t );     }
  
  explicit mat4( const T * t )
    { setValue ( t ); }

  mat4( const T & m00, const T & m01, const T & m02, const T & m03,
        const T & m10, const T & m11, const T & m12, const T & m13,  
        const T & m20, const T & m21, const T & m22, const T & m23,
        const T & m30, const T & m31, const T & m32, const T & m33 )  
    {
      setValue( m00, m01, m02, m03,
                m10, m11, m12, m13,  
                m20, m21, m22, m23,
                m30, m31, m32, m33 );
    }

  mat4( const mat4<T> & m ) 
    {
      _m[0] = m[0]; _m[4] = m[4]; _m[8 ] = m[8 ]; _m[12] = m[12];  
      _m[1] = m[1]; _m[5] = m[5]; _m[9 ] = m[9 ]; _m[13] = m[13];  
      _m[2] = m[2]; _m[6] = m[6]; _m[10] = m[10]; _m[14] = m[14];  
      _m[3] = m[3]; _m[7] = m[7]; _m[11] = m[11]; _m[15] = m[15];  
    }

  mat4( const vec4<T> & v0, const vec4<T> & v1, const vec4<T> & v2, const vec4<T> v3  ) 
    {
      _m[0] = v0[0]; _m[4] = v0[1]; _m[8]  = v0[2]; _m[12] = v0[3];  
      _m[1] = v1[0]; _m[5] = v1[1]; _m[9]  = v1[2]; _m[13] = v1[3];  
      _m[2] = v2[0]; _m[6] = v2[1]; _m[10] = v2[2]; _m[14] = v2[3];  
      _m[3] = v3[0]; _m[7] = v3[1]; _m[11] = v3[2]; _m[15] = v3[3];  
    }

  const T * getValue( T * t = NULL ) const
    // return 1D array observing OpenGL matrix element order, aka, column-major
    {  
      if (t!=NULL) { 
        t[0] = _m[0]; t[4] = _m[4]; t[8 ] = _m[8 ]; t[12] = _m[12];
        t[1] = _m[1]; t[5] = _m[5]; t[9 ] = _m[9 ]; t[13] = _m[13];
        t[2] = _m[2]; t[6] = _m[6]; t[10] = _m[10]; t[14] = _m[14];
        t[3] = _m[3]; t[7] = _m[7]; t[11] = _m[11]; t[15] = _m[15];
      }   
      return _m; 
    }

  void setValue( const T * t )
    {  
      _m[0] = t[0]; _m[4] = t[4]; _m[8 ] = t[8 ]; _m[12] = t[12];
      _m[1] = t[1]; _m[5] = t[5]; _m[9 ] = t[9 ]; _m[13] = t[13]; 
      _m[2] = t[2]; _m[6] = t[6]; _m[10] = t[10]; _m[14] = t[14];
      _m[3] = t[3]; _m[7] = t[7]; _m[11] = t[11]; _m[15] = t[15];
    }
  
  void getValue(  T & m11, T & m12, T & m13, T & m14,
                  T & m21, T & m22, T & m23, T & m24,  
                  T & m31, T & m32, T & m33, T & m34, 
                  T & m41, T & m42, T & m43, T & m44 ) const
    {
      m11 = _m[0]; m12 = _m[4];  m13 = _m[8 ]; m14 = _m[12];
      m21 = _m[1]; m22 = _m[5];  m23 = _m[9 ]; m24 = _m[13];
      m31 = _m[2]; m32 = _m[6];  m33 = _m[10]; m34 = _m[14];
      m41 = _m[3]; m42 = _m[7];  m43 = _m[11]; m44 = _m[15];
    }

  void setValue(  const T & m00, const T & m01, const T & m02, const T & m03,
                  const T & m10, const T & m11, const T & m12, const T & m13,  
                  const T & m20, const T & m21, const T & m22, const T & m23,
                  const T & m30, const T & m31, const T & m32, const T & m33 )  
    {
      _m[0] = m00; _m[4] = m01;  _m[8 ] = m02; _m[12] = m03; 
      _m[1] = m10; _m[5] = m11;  _m[9 ] = m12; _m[13] = m13;   
      _m[2] = m20; _m[6] = m21;  _m[10] = m22; _m[14] = m23; 
      _m[3] = m30; _m[7] = m31;  _m[11] = m32; _m[15] = m33; 
    }

  vec4<T> getRowVector( int i ) const
    { return vec4<T>(_m[i], _m[i+4], _m[i+8], _m[i+12]); }

  void setRowVector( int i, const vec4<T> & v )
    { _m[i] = v[0]; _m[i+4] = v[1]; _m[i+8] = v[2]; _m[i+12] = v[3]; }

  void getRowValue( int i , T * t ) const
    { t[0] = _m[i]; t[1] = _m[i+4]; t[2] = _m[i+8]; t[3] = _m[i+12]; }

  void setRowValue( int i , const T * t )
    { _m[i] = t[0]; _m[i+4] = t[1]; _m[i+8] = t[2]; _m[i+12] = t[3]; }
    
  vec4<T> getColVector( int j ) const
    { return vec4<T>(&_m[j<<2]); }

  void setColVector( int j, const vec4<T> & v )
    { _m[j<<2] = v[0]; _m[(j<<2)+1] = v[1]; _m[(j<<2)+2] = v[2]; _m[(j<<2)+3] = v[3]; }

  void getColValue( int j , T * t ) const
    { t[0] = _m[j<<2]; t[1] = _m[(j<<2)+1]; t[2] = _m[(j<<2)+2]; t[3] = _m[(j<<2)+3]; }

  void setColValue( int j, const T * t )
    { _m[j<<2] = t[0]; _m[(j<<2)+1] = t[1]; _m[(j<<2)+2] = t[2]; _m[(j<<2)+3] = t[3]; }

  void reset( const T & t )
  {
      _m[0] = t; _m[4] = t; _m[8 ] = t; _m[12]  = t;
      _m[1] = t; _m[5] = t; _m[9 ] = t; _m[13]  = t; 
      _m[2] = t; _m[6] = t; _m[10] = t; _m[14]  = t; 
      _m[3] = t; _m[7] = t; _m[11] = t; _m[15]  = t; 
  }

  void zero()
    {  reset( 0 ); }
  
  void negate()
    {
      _m[0] = -_m[0]; _m[4] = -_m[4]; _m[8 ] = -_m[8 ]; _m[12]  = -_m[12];  
      _m[1] = -_m[1]; _m[5] = -_m[5]; _m[9 ] = -_m[9 ]; _m[13]  = -_m[13];  
      _m[2] = -_m[2]; _m[6] = -_m[6]; _m[10] = -_m[10]; _m[14]  = -_m[14];  
      _m[3] = -_m[3]; _m[7] = -_m[7]; _m[11] = -_m[11]; _m[15]  = -_m[15];  
    }
  
  bool isZero() const
    {  
      return  is_zero(_m[0 ]) && is_zero(_m[1 ]) && is_zero(_m[2 ]) && is_zero(_m[3 ]) &&
              is_zero(_m[4 ]) && is_zero(_m[5 ]) && is_zero(_m[6 ]) && is_zero(_m[7 ]) &&
              is_zero(_m[8 ]) && is_zero(_m[9 ]) && is_zero(_m[10]) && is_zero(_m[11]) &&
              is_zero(_m[12]) && is_zero(_m[13]) && is_zero(_m[14]) && is_zero(_m[15]) ;
    }

  bool agree( const mat4<T> & m) const
    {  return (*this - m).isZero();  }

  void identity()
    { 
      _m[0] = T(1); _m[4] = 0;     _m[8 ] = 0;    _m[12] = 0;  
      _m[1] = 0;    _m[5] = T(1);  _m[9 ] = 0;    _m[13] = 0;  
      _m[2] = 0;    _m[6] = 0;     _m[10] = T(1); _m[14] = 0;  
      _m[3] = 0;    _m[7] = 0;     _m[11] = 0;    _m[15] = T(1); 
    }
  
  void transpose()
    { 
      T t;
      t = a01; a01 = a10; a10 = t;
      t = a12; a12 = a21; a21 = t;
      t = a20; a20 = a02; a02 = t;
      t = a03; a03 = a30; a30 = t;
      t = a13; a13 = a31; a31 = t;
      t = a23; a23 = a32; a32 = t;
    }
/*
  T minor(int r0, int r1, int r2,  int c0, int c1, int c2) const
    // calculate the minor(r,c), where r = {0,1,2,3}-{r0,r1,r1}, c = {0,1,2,3}-{c0,c1,c2}
    { 
#define R (*this)
      return R(r0,c0) * ( R(r1,c1)*R(r2, c2) - R(r1,c2)*R(r2, c1) )+
             R(r0,c1) * ( R(r1,c2)*R(r2, c0) - R(r1,c0)*R(r2, c2) )+
             R(r0,c2) * ( R(r1,c0)*R(r2, c1) - R(r1,c1)*R(r2, c0) );
#undef R
    }

  T determinant() const
    { 
      return + _m[0]*minor(1,2,3, 1,2,3) 
             - _m[1]*minor(0,2,3, 1,2,3) 
             + _m[2]*minor(0,1,3, 1,2,3)
             - _m[3]*minor(0,1,2, 1,2,3); 
    }

  void adjoint()
    {
      mat4<T> m(*this);
#define M(r0,r1,r2, c0,c1,c2)  m.minor(r0,r1,r2, c0,c1,c2)
      _m[0] =  M(1,2,3, 1,2,3); _m[4] = -M(0,2,3, 1,2,3); _m[8 ] =  M(0,1,3, 1,2,3); _m[12] = -M(0,1,2, 1,2,3);
      _m[1] = -M(1,2,3, 0,2,3); _m[5] =  M(0,2,3, 0,2,3); _m[9 ] = -M(0,1,3, 0,2,3); _m[13] =  M(0,1,2, 0,2,3);
      _m[2] =  M(1,2,3, 0,1,3); _m[6] = -M(0,2,3, 0,1,3); _m[10] =  M(0,1,3, 0,1,3); _m[14] = -M(0,1,2, 0,1,3);
      _m[3] = -M(1,2,3, 0,1,2); _m[7] =  M(0,2,3, 0,1,2); _m[11] = -M(0,1,3, 0,1,2); _m[15] =  M(0,1,2, 0,1,2);
#undef M
    }
*/
  bool inverse()
    {
     /*
	  T det = determinant();

      if (!is_zero(det)) {
        this->adjoint();
        *this *= T(1)/det;
        return true;
      }
      */
      return false;
    }

  bool isSingular () const
    {
      /*
	  T det = determinant();
      return is_zero( det );
      */
	  return false;
    }

  T rowDotCol(int r, mat4<T> rhs, int c) const
    { 
      return  (*this)(r,0) * rhs(0,c) + 
              (*this)(r,1) * rhs(1,c) + 
              (*this)(r,2) * rhs(2,c) + 
              (*this)(r,3) * rhs(3,c) ;
    }

  void print() const
    {  
      std::cout << "[ " << _m[0] << " " << _m[4] <<" " << _m[8 ] << " " << _m[12] << " ]\n"
                << "[ " << _m[1] << " " << _m[5] <<" " << _m[9 ] << " " << _m[13] << " ]\n"
                << "[ " << _m[2] << " " << _m[6] <<" " << _m[10] << " " << _m[14] << " ]\n" 
                << "[ " << _m[3] << " " << _m[7] <<" " << _m[11] << " " << _m[15] << " ]\n";
    }
  
  // operator functions
  
  T & operator () ( int r, int c )
    { return _m[(c<<2)+r]; }

  const T & operator () ( int r, int c ) const
    { return _m[(c<<2)+r]; }

  T & operator [] ( int i )
    { return _m[i]; }
  
  const T & operator [] ( int i ) const
    { return _m[i]; }

  mat4<T> & operator = ( const mat4<T> & m )
    { 
      _m[0] = m[0]; _m[4] = m[4]; _m[8 ] = m[8 ]; _m[12] = m[12];  
      _m[1] = m[1]; _m[5] = m[5]; _m[9 ] = m[9 ]; _m[13] = m[13];  
      _m[2] = m[2]; _m[6] = m[6]; _m[10] = m[10]; _m[14] = m[14];  
      _m[3] = m[3]; _m[7] = m[7]; _m[11] = m[11]; _m[15] = m[15];  
      return (*this); 
    }
  
  mat4<T> & operator += ( const mat4<T> & m )
    { 
      _m[0] += m[0]; _m[4] += m[4]; _m[8 ] += m[8 ];  _m[12] += m[12];  
      _m[1] += m[1]; _m[5] += m[5]; _m[9 ] += m[9 ];  _m[13] += m[13];  
      _m[2] += m[2]; _m[6] += m[6]; _m[10] += m[10]; _m[14] += m[14];  
      _m[3] += m[3]; _m[7] += m[7]; _m[11] += m[11]; _m[15] += m[15];  
      return (*this); 
    }
  
  mat4<T> & operator -= ( const mat4<T> & m )
    { 
      _m[0] -= m[0]; _m[4] -= m[4]; _m[8 ] -= m[8 ]; _m[12] -= m[12];  
      _m[1] -= m[1]; _m[5] -= m[5]; _m[9 ] -= m[9 ]; _m[13] -= m[13];  
      _m[2] -= m[2]; _m[6] -= m[6]; _m[10] -= m[10]; _m[14] -= m[14];  
      _m[3] -= m[3]; _m[7] -= m[7]; _m[11] -= m[11]; _m[15] -= m[15];  
      return (*this); 
    }

  mat4<T> & operator *= ( const mat4<T> & rhs )
    { 
      mat4<T> lhs(*this);
#define D(r, c)  lhs.rowDotCol(r, rhs, c)
      _m[0] = D(0, 0); _m[4] = D(0, 1); _m[8 ] = D(0, 2); _m[12] = D(0, 3); 
      _m[1] = D(1, 0); _m[5] = D(1, 1); _m[9 ] = D(1, 2); _m[13] = D(1, 3); 
      _m[2] = D(2, 0); _m[6] = D(2, 1); _m[10] = D(2, 2); _m[14] = D(2, 3); 
      _m[3] = D(3, 0); _m[7] = D(3, 1); _m[11] = D(3, 2); _m[15] = D(3, 3); 
#undef D
      return (*this); 
    }

  mat4<T> & operator *= ( T t )
    { 
      _m[0] *= t;  _m[4] *= t;  _m[8 ] *= t;  _m[12] *= t;  
      _m[1] *= t;  _m[5] *= t;  _m[9 ] *= t;  _m[13] *= t;  
      _m[2] *= t;  _m[6] *= t;  _m[10] *= t;  _m[14] *= t;  
      _m[3] *= t;  _m[7] *= t;  _m[11] *= t;  _m[15] *= t;  
      return (*this); 
    }
  
  mat4<T> & operator /= ( T t )
    { 
      if (t != T(0)) {
        T rt = T(1)/t;
        _m[0] *= rt;  _m[4] *= rt;  _m[8 ] *= rt;  _m[12] *= rt;  
        _m[1] *= rt;  _m[5] *= rt;  _m[9 ] *= rt;  _m[13] *= rt;  
        _m[2] *= rt;  _m[6] *= rt;  _m[10] *= rt;  _m[14] *= rt;  
        _m[3] *= rt;  _m[7] *= rt;  _m[11] *= rt;  _m[15] *= rt;  
      }
      return (*this); 
    }
  
  mat4<T> operator + () const
    { return mat4<T>(*this); }
  
  mat4<T> operator - () const
    { mat4<T> rm(*this); rm.negate(); return rm; }
  
  mat4<T> operator + ( const mat4<T> & m ) const
    { mat4<T> rm(*this); return rm += m; }
  
  mat4<T> operator - ( const mat4<T> & m ) const
    { mat4<T> rm(*this); return rm -= m; }

  mat4<T> operator * ( const mat4<T> & m ) const
    { mat4<T> rm(*this); return rm *= m; }

  vec4<T> operator * ( const vec4<T> & v ) const
    { 
      return vec4<T>( getRowVector(0).dot(v), 
                      getRowVector(1).dot(v), 
                      getRowVector(2).dot(v), 
                      getRowVector(3).dot(v) ); 
    }

  mat4<T> operator * ( T t ) const
    { mat4<T> rm(*this); return rm *= t; }
  
  mat4<T> operator / ( T t ) const
    { mat4<T> rm(*this); return rm /= t; }

  bool operator == ( const mat4<T> & m ) const
    {  
      return  _m[0]==m[0] && _m[4]==m[4] && _m[8 ]==m[8 ] && _m[12]==m[12] && 
              _m[1]==m[1] && _m[5]==m[5] && _m[9 ]==m[9 ] && _m[13]==m[13] && 
              _m[2]==m[2] && _m[6]==m[6] && _m[10]==m[10] && _m[14]==m[14] && 
              _m[3]==m[3] && _m[7]==m[7] && _m[11]==m[11] && _m[15]==m[15] ; 
    }

  bool operator != ( const mat4<T> & m) const
    {  return !(*this == m);  }
};

// mat4<T> operator functions
		
template <typename T> inline
  mat4<T> operator * ( T t, const mat4<T> & m )
{ 
  return m*t; 
}

template <typename T> inline
  vec4<T> operator * ( const vec4<T> & v, const mat4<T> & m )
{ 
  return vec4<T>( v.dot(m.getColVector(0)), 
                  v.dot(m.getColVector(1)), 
                  v.dot(m.getColVector(2)), 
                  v.dot(m.getColVector(3)) ); 
}

// mat4<T> template functions

template <typename T> inline
  mat4<T>  transpose( const mat4<T> & m )  
{
  mat4<T> rm(m);
  rm.transpose();
  return rm;
}

template <typename T> inline
  mat4<T>  inverse( const mat4<T> & m )  
{
  mat4<T> rm(m);
  rm.inverse();
  return rm;
}


////
// Export:
//   2D, 3D, and 4D vector types
//   3x3 and 4x4 matrix types
////

typedef vec2<float> Vector2f;
typedef vec3<float> Vector3f;
typedef vec4<float> Vector4f;
typedef mat3<float> Matrix3f;
typedef mat4<float> Matrix4f;

typedef vec2<double> Vector2d;
typedef vec3<double> Vector3d;
typedef vec4<double> Vector4d;
typedef mat3<double> Matrix3d;
typedef mat4<double> Matrix4d;

typedef vec2<float> Vector2;
typedef vec3<float> Vector3;
typedef vec4<float> Vector4;
typedef mat3<float> Matrix3;
typedef mat4<float> Matrix4;

// explicit export template classes 
// only float based are exported to save code length since we only use float here
//template class VSL_DLL_SPEC vec2<float>;
//template class VSL_DLL_SPEC vec3<float>;
//template class VSL_DLL_SPEC vec4<float>;
//template class VSL_DLL_SPEC mat3<float>;
//template class VSL_DLL_SPEC mat4<float>;

#pragma warning( pop )

}

#endif



