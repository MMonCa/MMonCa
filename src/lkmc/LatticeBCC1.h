/*
 * LatticeBCC1.h
 *
 *  Created on: Aug 17, 2012
 *      Author: ignacio.martin@imdea.org
 */

#ifndef LATTICEBCC1_H_
#define LATTICEBCC1_H_

#include "Lattice.h"

namespace LKMC {

class LatticeBCC1: public LKMC::Lattice {
public:
	LatticeBCC1(Kernel::Domain *p, const LatticeParam *, Kernel::M_TYPE mt);
	virtual ~LatticeBCC1();
	LatticeParam::TYPE getType() const { return LatticeParam::BCC; }
	void fill(const Kernel::Coordinates &mc, const Kernel::Coordinates &Mc,
			std::vector<LatticeInformation> &pos, bool onlyFrame) const;
};

} /* namespace LKMC */
#endif /* LATTICEBCC1_H_ */
