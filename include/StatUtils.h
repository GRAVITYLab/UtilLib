#ifndef _STATUTILS_H  
#define _STATUTILS_H  

#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cassert>

#include "vslLinear.h"

using namespace std;

namespace utillib
{
	template <class T>
	class StatUtils
	{
	private:

	public:
		static T _ulib_stutils_clamp( T val, T min, T max );
		static T _ulib_stutils_Min( T *data, int numelements );
		static T _ulib_stutils_Max( T *data, int numelements );
		static T* _ulib_stutils_MinMaxRange( T *data, int numelements );

		static T _ulib_stutils_sum( T* array, int len );
		static T _ulib_stutils_prod( T* array, int len );
		static void _ulib_stutils_fill( T* array, int len, T val );

		static void _ulib_stutils_addArrays( T* one, T* two, T* sum, int len );
		static void _ulib_stutils_subtractArrays( T* one, T* two, T* sub, int len );
		static void _ulib_stutils_addArrayScalar( T* array, T scalar, int len );
		static void _ulib_stutils_subtractArrayScalar( T* array, T scalar, int len );
		static void _ulib_stutils_mulArrayScalar( T* array, T scalar, int len );

		static float _ulib_stutils_Lerp( float a, float b, float f );

		static T _ulib_stutils_Mean( T *data, int numelements );
		static T _ulib_stutils_StdDev( T *data, int numelements );

		static T* _ulib_stutils_normalizeMinMax( T *data, int numelements );
		static void _ulib_stutils_normalizeMuSigma( T *data, int numelements );
	};

	template <class T>
	T StatUtils<T>::_ulib_stutils_clamp( T val, T min, T max )
	{
		if( val < min ) return min;
		if( val > max ) return max;
		return val;
	}

	template <class T>
	T StatUtils<T>::_ulib_stutils_Min( T *data, int numelements )
	{
		T min = *(data);
		for( int i=0; i<numelements; i++ )
		{
			if( *(data + i) < min ) min = *(data + i);
		}
		return min;
	}// end function

	template <class T>
	T StatUtils<T>::_ulib_stutils_Max( T *data, int numelements )
	{
		T max = *(data);
		for( int i=0; i<numelements; i++ )
		{
			if( *(data + i) > max ) max = *(data + i);
		}
		return max;

	}// end function

	template <class T>
	T StatUtils<T>::_ulib_stutils_sum( T* array, int len )
	{
		T sum = 0;
		for( int i=0; i<len; i++ )
			sum += array[i];
		return sum;

	}// end function

	template <class T>
	T StatUtils<T>::_ulib_stutils_prod( T* array, int len )
	{
		T product = 1;
		for( int i=0; i<len; i++ )
			product *= array[i];
		return product;

	}// end function

	template <class T>
	void StatUtils<T>::_ulib_stutils_fill( T* array, int len, T val )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = val;
		}
	}

	template <class T>
	void StatUtils<T>::_ulib_stutils_addArrays( T* one, T* two, T* sum, int len )
	{
		assert( sum != NULL );
		for( int i=0; i<len; i++ )
		{
			sum[i] = one[i]+two[i];
		}
	}

	template <class T>
	void StatUtils<T>::_ulib_stutils_subtractArrays( T* one, T* two, T* sub, int len )
	{
		assert( sub != NULL );
		for( int i=0; i<len; i++ )
		{
			sub[i] = one[i]-two[i];
		}
	}

	template <class T>
	void StatUtils<T>::_ulib_stutils_addArrayScalar( T* array, T scalar, int len )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = array[i]+scalar;
		}
	}

	template <class T>
	void StatUtils<T>::_ulib_stutils_subtractArrayScalar( T* array, T scalar, int len )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = array[i]-scalar;
		}
	}

	template <class T>
	void StatUtils<T>::_ulib_stutils_mulArrayScalar( T* array, T scalar, int len )
	{
		for( int i=0; i<len; i++ )
		{
			array[i] = array[i]*scalar;
		}
	}

	template <class T>
	T* StatUtils<T>::_ulib_stutils_MinMaxRange( T *data, int numelements )
	{
		T *ret = new T[3];
		ret[0] = ret[1] = *data;
		ret[2] = 0;

		for( int i=0; i<numelements; i++ )
		{
			if( *(data + i) < ret[0] ) ret[0] = *(data + i);
			if( *(data + i) > ret[1] ) ret[1] = *(data + i);
		}
		ret[2] = ret[1] - ret[0];
		
		return ret;

	}// end function
	
	template <class T>
	float StatUtils<T>::_ulib_stutils_Lerp( float a, float b, float f )
	{
		if( a<b )	return a + f * fabs( a-b );
		return a - f * fabs( a-b );
	}

	template <class T>
	T StatUtils<T>::_ulib_stutils_Mean( T *data, int numelements )
	{
		T sum = 0;
		for( int i=0; i<numelements; i++ )
			sum += *(data + i);

		return ( sum / (T)numelements );
	}// end function

	template <class T>
	T StatUtils<T>::_ulib_stutils_StdDev( T *data, int numelements )
	{
		T mean = StatUtils::_ulib_stutils_Mean( data, numelements );

		T score = 0;
		for( int i=0; i<numelements; i++ )
			score += (*(data + i) - mean) * (*(data + i) - mean);
		
		return sqrt( score / (T)numelements );
	}// end function

	template <class T>
	T* StatUtils<T>::_ulib_stutils_normalizeMinMax( T *data, int numelements )
	{
		T* out = StatUtils::_ulib_stutils_MinMaxRange( data, numelements );
		cout << "out here" << endl;
		cout << out[0] << "," << out[1] << "," << out[2] << endl;
		
		T* ret = new T[3];
		// Return all zeros, if all values are same
		if( out[2] == 0 )
		{
			cout << "here" << endl;
			for( int i=0; i<numelements; i++ )
				data[i] = 0;

			ret[0] = ret[1] = ret[2] = 0;
			return ret;
		}

		// Otherwise normalize each value of the array
		for( int i=0; i<numelements; i++ )
		{
			//if( data[i] != 0 ) cout << data[i] << endl;
			data[i] = (data[i] - out[0]) / (T)out[2];

			if( i == 0 )	ret[0] = ret[1] = data[i];
			else
			{
				if( data[i] < ret[0] )	ret[0] = data[i];
				if( data[i] > ret[1] )	ret[1] = data[i];
			}
		}

		return ret;

	}// end function

	template <class T>
	void StatUtils<T>::_ulib_stutils_normalizeMuSigma( T *data, int numelements )
	{
		T mean = StatUtils::_ulib_stutils_Mean( data, numelements );
		T sigma = StatUtils::_ulib_stutils_StdDev( data, numelements );

		for( int i=0; i<numelements; i++ )
			*(data + i) = (*(data + i) - mean) / sigma;

	}// end function

}// end namespace

#endif
