/*
 * EpiGasParam.h
 *
 *  Created on: Oct 22, 2014
 *      Author: ignacio.martin@imdea.org
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

#ifndef EPIGASPARAM_H_
#define EPIGASPARAM_H_

#include <map>
#include <string>
#include <tcl.h>
#include "kernel/ParticleType.h"
#include "kernel/Material.h"
#include "kernel/Domain.h"

namespace LKMC {

const unsigned MAX_EPI_GASES = 3;  //                    +1                 +2                  +3
enum LkmcEvent { LE_FINAL_MIGRATES = MAX_EPI_GASES, LE_FINAL_REMOVED, LE_PRECUR_TO_FINAL, LE_PRECUR_REMOVED, MAX_LKMC_EVENT };

struct EpiGasParam
{
	EpiGasParam();

	float _partialPress[Kernel::MAX_IMPURITIES];
	float _prefPrecursor[Kernel::MAX_IMPURITIES]; //here to be updated with temperature and partial pressures
	Kernel::P_TYPE _toPType[MAX_EPI_GASES];
	void processLine(Tcl_Interp *pTcl, const Kernel::Domain *, Kernel::M_TYPE, double T, const std::map<std::string, float> &pressures);

private:
	void initialize();
};

} /* namespace LKMC */

#endif /* EPIGASPARAM_H_ */
