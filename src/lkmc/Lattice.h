/*
 *
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

#ifndef LKMCLATTICE_H
#define LKMCLATTICE_H

#include "kernel/Coordinates.h"
#include "LatticeParam.h"
#include <vector>
#include <set>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace Kernel
{
	class Domain;
}

namespace LKMC 
{	
class LatticeAtom;

class Lattice
{
public:
	struct LatticeInformation
	{
		Kernel::Coordinates _coords;
		unsigned            _type;
		int                 _orientation;
	};

	Lattice(Kernel::Domain *p, const LatticeParam *, Kernel::M_TYPE mt);
    virtual ~Lattice() {}

    virtual void fill(const Kernel::Coordinates &m, const Kernel::Coordinates &M,
              std::vector<LatticeInformation> &, bool onlyFrame) const = 0;
    virtual LatticeParam::TYPE getType() const = 0;

    Kernel::Coordinates toRotated  (const Kernel::Coordinates &o)  const;
    Kernel::Coordinates fromRotated(const Kernel::Coordinates &n) const;

    void getMatR(ublas::matrix<double> &m) const;
    void getInvR(ublas::matrix<double> &m) const;

protected:
    Kernel::Domain *_pDomain;
    float _latticeParameter[2];  //[alloy]
    double _matR[3][3];
    double _invR[3][3];
    Kernel::M_TYPE _mt;
    
    void setMatrices(float Mi, float Mj, float Mk, float Fi, float Fj, float Fk);
    void getIndices(const Kernel::Coordinates &mc, const Kernel::Coordinates &Mc, int m[3], int M[3]) const;

    void toLatticeUnits(int &mi, int &mj, int &mz, const Kernel::Coordinates &) const;
    //If the atom specified by c fills in mc and Mc, then it is filled in pos
    void tryInBox(std::vector<LatticeInformation> &pos, const Kernel::Coordinates &mc, const Kernel::Coordinates &Mc,
    	int orientation, unsigned type, const Kernel::Coordinates &c) const;
};

Lattice * readLattice(Kernel::Domain *pD, const LatticeParam *p, Kernel::M_TYPE mt);

}

#endif
