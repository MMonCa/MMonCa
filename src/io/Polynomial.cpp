/*
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

/* This structure contains a set of regions, each with boundaries, and a polynomia
 * between them. The polynomia are defined as a bunch of N numbers (lets say 3,2,...,4)
 * defining 3**N-1 + 2**(N-2) + ... + 4
 */

#include "Polynomial.h"
#include <math.h>

namespace IO {

Polynomial::Polynomial(std::vector<std::vector <double> > vPoly, std::vector <double> bounds,
		std::vector<double> centered)
{
	_precission = 1001;
	_derivated = false;
	_vPoly = vPoly;
	_bounds = bounds;
	_centered = centered;
	derive();
	double x;
	for(unsigned idx = 0; idx <= _precission; idx++)
	{
		x = double(idx) / double(_precission);
		_px.push_back(evaluate(_vPoly, x));
		_dpx.push_back(evaluate(_vdPoly, x));
	}
}

Polynomial::Polynomial()
{
	_bounds.push_back(0);
	_centered = _bounds;
	_vPoly.push_back(_bounds);
	_vdPoly.push_back(_bounds);
	_precission = 1001;
	_derivated = true;
}

double Polynomial::getValue(double x)
{
	return evaluate(_vPoly, x);
}

double Polynomial::getFirstDerValue(double x)
{
	return evaluate(_vdPoly, x);
}

double Polynomial::getIntValue(double x)
{
	return interpolate(_px, x);
}

double Polynomial::getIntFirstDerValue(double x)
{
	return interpolate(_dpx, x);
}

double Polynomial::evaluate(const std::vector < std::vector<double> >& vV , double x)
{
	double value = 0;
	if (vV.size() == 0)
		return value;
		
	for(unsigned bnd = 0; bnd < (_bounds.size() - 1); bnd++)
		if(_bounds[bnd] <= x && x < _bounds[bnd + 1]) // I am on interval "bnd"
			for(unsigned grade = 0; grade < vV[bnd].size(); grade++)
				value += vV[bnd][grade] * pow(x - _centered[bnd], vV[bnd].size() - 1 - grade);
	return value;
}

double Polynomial::interpolate(const std::vector<double>& pts, double x)
{
        if(x == 0)
            return pts[0];
        if(x == 1)
            return pts[_precission];
	int idx = int(x * (_precission - 1));
	double xi, xf;
	xi = double(idx) / double(_precission - 1);
	xf = double(idx + 1) / double(_precission - 1);
	return pts[idx] + (x - xi) * (pts[idx + 1] - pts[idx]) / (xf - xi);
}

void Polynomial::derive()
{
	if(_derivated)
		return;
	for(const auto &pol : _vPoly)
	{
		std::vector <double> dervPol;
		if(pol.size())
			for (unsigned i = 0; i < (pol.size() - 1); i++)
			{
				int exp = (pol.size() - 1 - i);
				dervPol.push_back(exp * pol[i]);
			}
		else
			dervPol.push_back(0.0);
		_vdPoly.push_back(dervPol);
	}
	_derivated = true;
}

}
