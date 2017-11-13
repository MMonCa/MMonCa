/*
 * ArrheniusAlloys.h
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

#ifndef ARRHENIUSALLOYS_H_
#define ARRHENIUSALLOYS_H_

#include <vector>
#include "Arrhenius.h"

namespace Kernel{
class MeshElement;
}

namespace IO
{
struct ArrheniusAlloys {
	static const unsigned NVALUES = 10;

	Arrhenius _v[NVALUES+1];

	ArrheniusAlloys() {}
	ArrheniusAlloys(double p, double e);
	ArrheniusAlloys(const std::vector<float> &p, const std::vector<float> &e, const std::vector<float> &c);

	Arrhenius 		 operator()	(double) const;
	Arrhenius 		 operator()	(const Kernel::MeshElement *pME) const;
	ArrheniusAlloys& operator=  (const ArrheniusAlloys &);

	//well defined algebraic operations for pref*exp(-E/kT) at the same T
	ArrheniusAlloys  operator*  (const ArrheniusAlloys &a) const { ArrheniusAlloys tmp(*this); tmp*=a;	return tmp; }
	ArrheniusAlloys  operator/  (const ArrheniusAlloys &a) const { ArrheniusAlloys tmp(*this); tmp/=a;	return tmp; }
	ArrheniusAlloys& operator*= (const ArrheniusAlloys &);
	ArrheniusAlloys& operator/= (const ArrheniusAlloys &);

	ArrheniusAlloys  operator*	(double a) const { ArrheniusAlloys tmp(*this); tmp*=a;	return tmp; }
	ArrheniusAlloys  operator/	(double a) const { ArrheniusAlloys tmp(*this); tmp/=a;	return tmp; }
	ArrheniusAlloys& operator*=	(double);
	ArrheniusAlloys& operator/= (double);

	bool             canDivide() const;
	void 			 addEner	(ArrheniusAlloys&);  //FIXME: It might not have physical sense. Use *= instead.
	void 			 addEner	(double); //to add barriers

private:
	Arrhenius 		 interpolate(const std::vector <float> &p, const std::vector <float> &e, const std::vector <float> &c, float comp) const;
};

inline std::ostream & operator << (std::ostream &is, const ArrheniusAlloys &arr)
{
	for (unsigned i = 0; i <= arr.NVALUES; i++)
		is << float(i) / arr.NVALUES << "% " << arr._v[i]._pref << "*exp(" << -arr._v[i]._ener << "/kBT)\n";
	return is;
}

}
#endif /* ARRHENIUSALLOYS_H_ */
