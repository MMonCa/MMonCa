/*
 * Arrhenius.h
 *
 * Created on: Mar 18, 2011
 * Author: ignacio.martin@imdea.org
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

#ifndef ARRHENIUS_H_
#define ARRHENIUS_H_

#include <ostream>
#include <cmath>
#include <vector>           // Include for Arrhenius(vector, vector)

namespace IO
{
struct Arrhenius
{
	double _pref;
	double _ener;

	double getRate(double kT) const	{return _pref * std::exp(-_ener/kT);}
	Arrhenius() : _pref(0), _ener(5) { }
	Arrhenius(double p, double e) : _pref(p), _ener(e) { }

	Arrhenius   operator* (const Arrhenius &arg) const { Arrhenius a(*this);  a *= arg; return a; }
	Arrhenius   operator/ (const Arrhenius &arg) const { Arrhenius a(*this);  a /= arg; return a; }
	Arrhenius & operator*=(const Arrhenius &arg) { _pref *= arg._pref; _ener += arg._ener; return *this; }
	Arrhenius & operator/=(const Arrhenius &arg) { _pref /= arg._pref; _ener -= arg._ener; return *this; }

	Arrhenius   operator* (double p) const { Arrhenius a(*this);  a._pref *= p; return a; }
	Arrhenius   operator/ (double p) const { Arrhenius a(*this);  a._pref /= p; return a; }
	Arrhenius & operator*=(double p) { _pref *= p; return *this; }
	Arrhenius & operator/=(double p) { _pref /= p; return *this; }
};


inline std::ostream & operator << (std::ostream &is, const Arrhenius &a)
{
	is << a._pref << "*exp(-" << a._ener << "/kBT)";
	return is;
}

}


#endif /* ARRHENIUS_H_ */
