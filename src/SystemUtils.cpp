/*
 * SystemUtils.cpp
 *
 *  Created on: Aug 7, 2013
 *      Author: abon
 */

#include "SystemUtils.h"

namespace utillib
{

	/**
	 * Start Time computation.
	 */
	double
	SystemUtils::startTimer()
	{
		return MPI_Wtime();
	}// end function

	/**
	 * Complete Time computation.
	 * @param startTime Recored time at start point
	 * @return computed time in seconds from start to end.
	 */
	double
	SystemUtils::endTimer( double startTime )
	{
		double endTime = MPI_Wtime();
		return ( endTime - startTime );

	}// end function

}

