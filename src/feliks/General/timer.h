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
/*
 * timer.h
 *
 * i. romero, october 2002
 * revised sept 2003
 * converted to c++ march 2006
 *
 * in this module a sort of stopwatch is implemented, in such a
 * way that many of them can be created and measure different 
 * times. Time is measured in miliseconds.
 *
 * The regular use should be something like
 *
 *    timer   w;    			// the declaration
 *    double  sec;  			// elapsed seconds
 *    double  msec;				// the elapsed miliseconds;
 *        ...
 *    w.start();				// start the clock
 *        xxx
 *        xxx
 *    w.stop();
 *    sec =  w.seconds();		// obtain the time from the starting point
 *    msec = w.miliseconds();
 *    w.start();
 *    w.stop();
 *    sec =  w.seconds();		// obtain the time from the starting point
 *    w.reset();
 *    w.start();
 *        ...
 *        ...
 */


#ifndef _timer_h
#define _timer_h

#include <sys/time.h>

class timer
{
private:
	long    oldclocks;
	clock_t startTime;
	clock_t stopTime;
	bool    ticking;
	struct timeval	startWallTime;
	struct timeval  stopWallTime;
	long   int  oldWallTimeMsecs;
	
	static long ctimeElapsedMsecs(const struct timeval& tstart, const struct timeval& tend);

public:

	timer();
	double  miliseconds();
	void    start();
	void    stop();
	void    reset();
	double  seconds();
	long	wallTimeMiliSeconds() const;
};

#endif
