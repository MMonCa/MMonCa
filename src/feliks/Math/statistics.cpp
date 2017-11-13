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
#include "Math/statistics.h"

#include <cmath>
#include <ctime>
#include <cstdlib>


bool random::initialized = false;


int  random :: discreteUniform(const int low, const int up)
{
    return ceil( uniform( (double) low, (double) up) );
}


double random :: uniform(const double low, const double up)
{
	if (!initialized) random::initialize();
	int     random_int = rand();
	double  random_dou = ( static_cast<double_t>(random_int))/ static_cast<double_t>(RAND_MAX);

	return random_dou * (up-low) + low;
}


double random :: weibull(const double shape, const double scale)
{
	if (!initialized) random::initialize();
	double m = - log(( (double) rand()/ (double) RAND_MAX) * (1.0-0.0) );
	return (scale * pow(m, 1.0/shape));
}


double random :: exponential(const double lambda)
{
	if (!initialized) random::initialize();
	double m = - log(1.0  - ( (double) rand()/ (double) RAND_MAX) );
	return (m/lambda);
}



void random :: initialize()
{
	srand((unsigned)time(0)); 
	initialized = true;
}
