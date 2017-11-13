/*
 * LatticeFCC1.h
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

#ifndef LATTICEFCC1_H_
#define LATTICEFCC1_H_

#include "Lattice.h"

namespace LKMC {

class LatticeFCC1: public LKMC::Lattice {
public:
	LatticeFCC1(Kernel::Domain *p, const LatticeParam *, Kernel::M_TYPE mt);
	virtual ~LatticeFCC1();
	LatticeParam::TYPE getType() const { return LatticeParam::FCC; }
	void fill(const Kernel::Coordinates &mc, const Kernel::Coordinates &Mc,
				std::vector<LatticeInformation> &pos, bool onlyFrame) const;
};

} /* namespace LKMC */
#endif /* LATTICEFCC1_H_ */
