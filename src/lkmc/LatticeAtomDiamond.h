/*
 * LatticeAtomDiamond.h
 *
 *  Created on: Oct 10, 2011
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

#ifndef LATTICEATOMDIAMOND_H_
#define LATTICEATOMDIAMOND_H_

#include "LatticeAtom.h"
#include "LKMCModel.h"
#include "LKMCMode.h"
#include <vector>
#include <map>

using namespace boost::numeric;

namespace LKMC {

class LatticeAtomDiamond: public LatticeAtom
{
public:
	LatticeAtomDiamond(LASTATE, Kernel::Domain *p, Kernel::M_TYPE, Kernel::P_TYPE, const Kernel::Coordinates &c, unsigned orientation);
	LatticeAtomDiamond(std::istream &);
	virtual ~LatticeAtomDiamond() {}

	virtual LatticeParam::TYPE getLatticeType() const { return LatticeParam::DIAMOND; }

	virtual void perform (Kernel::SubDomain *, unsigned eventType) = 0;
	virtual float getRate(unsigned eventType, float kT) const = 0;
	virtual void setIndex(unsigned eventType, int idx) = 0;
	virtual int  getIndex(unsigned eventType) const = 0;
	virtual unsigned getNEvents() const = 0;

	virtual void  setNeiDist(int i, float d2) { /*_neiDist[i] = std::sqrt(d2);*/ };
	virtual float getNeiDist(int i) const { return 0;/*_neiDist[i];*/ }
	virtual bool isDefective() const;
	unsigned getOrientation() const { return _orientation; }

	void restart(std::ostream &) const;
	static void toBeUpdated(LatticeAtomDiamond *, std::map<int, LatticeAtomDiamond *> &);

protected:
	template <typename AtomType>
	void createFromNeiList(Kernel::SubDomain *, LASTATE, const std::vector<std::pair<Kernel::Coordinates, unsigned> > &neis,
	    AtomType *pAmo, LKMC::LKMCMode, std::map<int, LatticeAtomDiamond *> &toUpdate) const;

	bool isAllowed(const Kernel::Coordinates &, const Kernel::MeshElement *, LKMCMode) const;
	unsigned _orientation;
};

}

#endif /* LATTICEATOMSI_H_ */
