/*
 * Polynomial.h
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain
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

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <vector>

namespace IO {

class Polynomial {
public:
	Polynomial(std::vector<std::vector <double> >, std::vector<double>, std::vector<double>);
	Polynomial();        
        
	double getIntValue(double);              // Interpolated Value (Fast)
	double getIntFirstDerValue(double);      // Interpolated Value (Fast)
	double getValue(double);              // Evaluated Value    (Slow)
	double getFirstDerValue(double);      // Evaluated Value    (Slow)

private:
	double evaluate(const std::vector< std::vector <double> > &, double);
	double interpolate(const std::vector<double> &, double);
	void  derive();
	std::vector< std::vector <double> > _vPoly;
	std::vector< std::vector <double> > _vdPoly;
	std::vector<double>                 _bounds;
	std::vector<double>                 _centered;
	bool							    _derivated;
	std::vector<double>				    _px;
	std::vector<double>                 _dpx;
	unsigned                            _precission; // Number of precalculated points (size of _px)
};

}
#endif /* POLYNOMIAL_H_ */

