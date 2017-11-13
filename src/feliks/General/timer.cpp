/* Finite Element Method Module
 *
 * Author: ignacio.romero@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain,
 *      and       Technical University of Madrid (UPM), Madrid, Spain
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */ 
/* in this module we implement a timer, an "object" that can act as 
   stopwatch, in a way that several of them can be running at the
   same time. 

   i. romero
   creation date: 14 novemeber 2002
modified: 16 september 2003
*/

#include <cstdlib>
#include <time.h>
#include "General/timer.h"



timer :: timer() :
	oldclocks(0),
	startTime(0),
	stopTime(0),
	ticking(false),
	oldWallTimeMsecs(0)
{
	startWallTime.tv_sec = startWallTime.tv_usec = 0;
	stopWallTime.tv_sec  = stopWallTime.tv_usec  = 0;
}


void timer :: start()
{
	oldclocks        += stopTime - startTime;
	oldWallTimeMsecs += ctimeElapsedMsecs(startWallTime, stopWallTime);
 
	startTime = clock();
	gettimeofday(&startWallTime, NULL);
	ticking   = true;
}

void timer :: stop()
{
	stopTime     = clock();
	gettimeofday(&stopWallTime, NULL);
	ticking      = false;
}


double timer :: miliseconds()
{
	double ms;
	
	if (ticking) ms =  (1000.0*(clock()  - startTime + oldclocks)) / ((double) CLOCKS_PER_SEC);
	else         ms =  (1000.0*(stopTime - startTime + oldclocks)) / ((double) CLOCKS_PER_SEC);
	
	return ms;
}


double timer :: seconds()
{
	double s;
	
	if (ticking) s =  (clock()  - startTime + oldclocks) / ((double) CLOCKS_PER_SEC);
	else         s =  (stopTime - startTime + oldclocks) / ((double) CLOCKS_PER_SEC);
	
	return s;
}


void timer :: reset()
{
	oldclocks = 0;
	oldWallTimeMsecs = 0;
	startTime = 0;
	stopTime  = 0;
	ticking   = false;
	startWallTime.tv_sec = startWallTime.tv_usec = 0;
	stopWallTime.tv_sec  = stopWallTime.tv_usec  = 0;
}


long timer :: wallTimeMiliSeconds() const
{
	long msecs;
	if (ticking) 
	{
		struct timeval tmp;
		gettimeofday(&tmp, NULL);
		msecs = ctimeElapsedMsecs(startWallTime, tmp) + oldWallTimeMsecs;
	}

	else
		msecs = ctimeElapsedMsecs(startWallTime, stopWallTime) + oldWallTimeMsecs;
	return msecs;
}




long timer :: ctimeElapsedMsecs(const struct timeval& tstart, const struct timeval& tend)
{
	return  (tend.tv_sec  - tstart.tv_sec) * 1000 + (tend.tv_usec - tstart.tv_usec) / 1000;
}

