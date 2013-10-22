/*
 * InterpUtils.cpp
 *
 *  Created on: Nov 7, 2012
 *      Author: abon
 */

#include "InterpUtils.h"

double
InterpUtils::trilinterpScalar( double* scalarArray, double* f )
{
	double t[2];
	double t_out;

	t[0] = bilinterpScalar( scalarArray, &f[0] );
	t[1] = bilinterpScalar( scalarArray + 4, &f[0] );
	t_out = linterpScalar( t, f[2] );

	return t_out;

}// end function

double
InterpUtils::bilinterpScalar( double* scalarArray, double* f )
{
	double t[2];
	double t_out;

	t[0] = linterpScalar( scalarArray, f[0] );
	t[1] = linterpScalar( scalarArray + 2, f[0] );
	t_out = linterpScalar( t, f[1] );

	return t_out;

}// end function

double
InterpUtils::linterpScalar( double* scalarArray, double f )
{
	double t_out = scalarArray[0] * f + scalarArray[1] * (1-f);
	return t_out;

}// end function

VECTOR3
InterpUtils::trilinterpVector( VECTOR3* vecArray, double* f )
{
	VECTOR3 t[2];
	VECTOR3 t_out;

	t[0] = bilinterpVector( vecArray, &f[0] );
	t[1] = bilinterpVector( vecArray + 4, &f[0] );
	t_out = linterpVector( t, f[2] );

	return t_out;

}// end function

VECTOR3
InterpUtils::bilinterpVector( VECTOR3* vecArray, double* f )
{
	VECTOR3 t[2];
	VECTOR3 t_out;

	t[0] = linterpVector( vecArray, f[0] );
	t[1] = linterpVector( vecArray + 2, f[0] );
	t_out = linterpVector( t, f[1] );

	return t_out;

}// end function

VECTOR3
InterpUtils::linterpVector( VECTOR3* vecArray, double f )
{
	VECTOR3 t_out( 0.0, 0.0, 0.0 );

	for( int i=0; i<3; i++ )
		t_out[i] = vecArray[0][i] * f + vecArray[1][i] * (1-f);

	return t_out;

}// end function

int
InterpUtils::index3DTo1DScalar( int x, int y, int z, int* size )
{
	return ( size[0] * size[1] * z + size[0] * y + x );

}// end function

int
InterpUtils::index3DTo1DVector( int x, int y, int z, int* size )
{
	return ( size[0] * size[1] * z + size[0] * y + x ) * 3;

}// end function



