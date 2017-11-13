/*
 * LatticeBCC1Param.h
 *
 *  Created on: Aug 17, 2012
 *      Author: ignacio.martin@imdea.org
 */

#ifndef LATTICEBCC1PARAM_H_
#define LATTICEBCC1PARAM_H_

#include "LatticeParam.h"
#include "io/Arrhenius.h"

namespace LKMC {

class LatticeBCC1Param: public LKMC::LatticeParam {
public:

	static const unsigned FIRST = 8;
	static const unsigned SECOND = 8;

	LatticeBCC1Param(const IO::ParameterManager *, const IO::FileParameters *, Kernel::M_TYPE);
	virtual ~LatticeBCC1Param();

	IO::Arrhenius _activation[FIRST+1][SECOND+1];
};

} /* namespace Domains */
#endif /* LATTICEBCC1PARAM_H_ */
