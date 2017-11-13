/*
 * ArrheniusAlloys.cpp
 *
 * Created on: Dec 20, 2012
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

#include "ArrheniusAlloys.h"
#include "io/Diagnostic.h"
#include "kernel/MeshElement.h"
#include "io/ParameterManager.h"

using std::vector;

namespace IO
{

ArrheniusAlloys::ArrheniusAlloys(double p, double e)
{
	for(unsigned i=0; i <= NVALUES; ++i)
	{
		_v[i]._pref = p;
		_v[i]._ener = e;
	}
}

ArrheniusAlloys::ArrheniusAlloys(const vector <float> &p, const vector <float> &e, const vector <float> &c)
{
	for(unsigned i=0; i <= NVALUES; ++i)
		_v[i] = interpolate(p, e, c, float(i) / NVALUES);
}

Arrhenius ArrheniusAlloys::interpolate(
		const vector <float> &vpref, const vector <float> &vener, const vector <float> &vx,
		float x) const
{
	Arrhenius arr;
	if (x == 0)
	{
		arr._pref = vpref[0];
		arr._ener = vener[0];
		return arr;
	}
        if (x == 1)
        {
                arr._pref = vpref.back();
		arr._ener = vener.back();
		return arr;
        }
        else
        {
                for(unsigned i = 0; i < vx.size(); i++)
                        if(x >= vx[i] && x <= vx[i + 1])
                        {
                                arr._ener = vener[i] + (vener[i + 1] - vener[i]) * (x - vx[i]) / (vx[i + 1] - vx[i]);  // Linear interpolation for the energy
                                arr._pref = pow(vpref[i], (vx[i + 1] - x) / (vx[i + 1] - vx[i])) *
                                        pow(vpref[i + 1], (x - vx[i]) / (vx[i + 1] - vx[i]));  // Algebraic interpolation for the prefactor                                
                        }
                return arr;
        }
}

ArrheniusAlloys& ArrheniusAlloys::operator=(const ArrheniusAlloys& arg)
{
  for (unsigned i = 0; i <= NVALUES; i++)
	  _v[i] = arg._v[i];
  return *this;
}


//energies added, prefactors multiplied
ArrheniusAlloys & ArrheniusAlloys::operator*= (const ArrheniusAlloys &arg)
{
	for (unsigned i = 0; i <= NVALUES; i++)
		_v[i] *= arg._v[i];
	return *this;
}

ArrheniusAlloys & ArrheniusAlloys::operator/= (const ArrheniusAlloys &arg)
{
	for (unsigned i = 0; i <= NVALUES; i++)
		_v[i] /= arg._v[i];
	return *this;
}

ArrheniusAlloys & ArrheniusAlloys::operator*=(double factor)
{
	for (unsigned i = 0; i <= NVALUES; i++)
		_v[i] *= factor;
	return *this;
}

ArrheniusAlloys & ArrheniusAlloys::operator/=(double factor)
{
	for (unsigned i = 0; i <= NVALUES; i++)
		_v[i] /= factor;
	return *this;
}

void ArrheniusAlloys::addEner(ArrheniusAlloys& arg)
{
	for(unsigned i = 0; i <= NVALUES; i++)
		_v[i]._ener += arg._v[i]._ener;
}

void ArrheniusAlloys::addEner(double ener)
{
	for (unsigned i = 0; i <= NVALUES; i++)
		_v[i]._ener += ener;
}

Arrhenius ArrheniusAlloys::operator()(const Kernel::MeshElement *pME) const
{
	return operator ()(pME->getBAtoms() != 0 ? pME->getEffectiveAlloyFraction() : 0);
}

Arrhenius ArrheniusAlloys::operator()(double x) const
{        
	if(x == 0)
		return _v[0];
	if(x == 1)
		return _v[NVALUES];

	unsigned im = floor(x * NVALUES);
	unsigned iM = ceil (x * NVALUES);
	assert(x <= 1.);
	float xm = float(im) / NVALUES;
	float diffx = 1./ NVALUES;

	return Arrhenius(
			_v[im]._pref + (_v[iM]._pref - _v[im]._pref) * (x - xm) / diffx ,
			_v[im]._ener + (_v[iM]._ener - _v[im]._ener) * (x - xm) / diffx );
}

bool ArrheniusAlloys::canDivide() const
{
	for (unsigned i = 0; i <= NVALUES; i++)
		if(_v[i]._pref == 0)
			return false;
	return true;
}

}
