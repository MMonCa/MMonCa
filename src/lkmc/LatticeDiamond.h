/*
 * LatticeDiamond.h
 *
 *  Created on: Aug 17, 2012
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

#ifndef LATTICEDIAMOND_H_
#define LATTICEDIAMOND_H_

#include "Lattice.h"

namespace Kernel { class SubDomain; }

namespace LKMC {

class LatticeDiamond: public LKMC::Lattice {
public:
	LatticeDiamond(Kernel::Domain *p, const LatticeParam *, Kernel::M_TYPE mt);
	virtual ~LatticeDiamond();

	virtual LatticeParam::TYPE getType() const  { return LatticeParam::DIAMOND; }
	virtual void fill(const Kernel::Coordinates &m, const Kernel::Coordinates &M,
	              std::vector<LatticeInformation> &, bool onlyFrame) const;
	unsigned getNeis2(Kernel::SubDomain *, int alloy, const Kernel::Coordinates &one, const Kernel::Coordinates &two,
	    		std::vector<std::pair<Kernel::Coordinates, unsigned> > &all) const;
	void getNeis1(Kernel::SubDomain *, int alloy, const Kernel::Coordinates &one, unsigned ori,
				std::vector<std::pair<Kernel::Coordinates, unsigned> > &normal,
				std::vector<std::pair<Kernel::Coordinates, unsigned> > &twins) const;
	unsigned getNeis1(Kernel::SubDomain *, int alloy, const Kernel::Coordinates &one, const Kernel::Coordinates &two,
	    		std::vector<std::pair<Kernel::Coordinates, unsigned> > &all,
	    		unsigned orient1, unsigned orient2) const;
private:
	 void setSystems();
	 float _twinProb;
	 float _distortion;
	 std::vector<std::vector<Kernel::Coordinates> > _systems[2];  //normal and alloy

};

class LatticeDiamond2: public LatticeDiamond
{
public:
	LatticeDiamond2(Kernel::Domain *p, const LatticeParam *lp, Kernel::M_TYPE mt) : LatticeDiamond(p, lp, mt) {}
	virtual LatticeParam::TYPE getType() const  { return LatticeParam::DIAMOND2; }
};

} /* namespace LKMC */
#endif /* LATTICEDIAMOND_H_ */
