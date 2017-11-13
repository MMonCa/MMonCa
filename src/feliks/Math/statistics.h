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
#ifndef _statistics_h
#define _statistics_h

class random
{

public:

	// uniform distribution
	static double uniform(const double low, const double up);
	static int    discreteUniform(const int low, const int up);
    
	// usage random::weibull(shape, scale);
	static double weibull(const double shape, const double scale);

	// usage random::exponential(lambda);
	static double exponential(const double lambda);

	
private:

	static bool initialized;
	static void	initialize();
};


#endif
