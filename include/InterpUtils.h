/*
 * InterpUtils.h
 *
 *  Created on: Nov 7, 2012
 *      Author: abon
 */

#ifndef INTERPUTILS_H_
#define INTERPUTILS_H_

#include "VectorMatrix.h"

class InterpUtils
{
private:

public:

	static double trilinterpScalar( double* scalarArray, double* f );
	static double bilinterpScalar( double* scalarArray, double* f );
	static double linterpScalar( double* scalarArray, double f );

	static VECTOR3 trilinterpVector( VECTOR3* vecArray, double* f );
	static VECTOR3 bilinterpVector( VECTOR3* vecArray, double* f );
	static VECTOR3 linterpVector( VECTOR3* vecArray, double f );

	static int index3DTo1DScalar( int x, int y, int z, int* size );
	static int index3DTo1DVector( int x, int y, int z, int* size );

};
#endif
/* INTERPUTILS_H_ */
