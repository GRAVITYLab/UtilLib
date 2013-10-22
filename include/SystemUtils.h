/*
 * SystemUtils.h
 *
 *  Created on: Aug 7, 2013
 *      Author: abon
 */

#ifndef SYSTEMUTILS_H_
#define SYSTEMUTILS_H_

#include "mpi.h"

using namespace std;

namespace utillib
{
	class SystemUtils
	{
	public:
		static double startTimer();
		static double endTimer( double startTime );
	};
}


#endif /* SYSTEMUTILS_H_ */
