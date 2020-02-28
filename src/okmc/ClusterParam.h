/*
 * ClusterParam.h
 *
 *  Created on: May 30, 2011
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

#ifndef CLUSTERPARAM_H_
#define CLUSTERPARAM_H_

#include "kernel/ParticleType.h"
#include "kernel/Material.h"
#include "io/ArrheniusAlloys.h"
#include "EDType.h"
#include <vector>
#include <set>
#include <map>
#include <string>

namespace IO { class ParameterManager; class FileParameters; }
namespace Kernel { class Domain; class MeshElement; class SubDomain; }
struct Tcl_Interp;

namespace OKMC {

class ClusterParam {
public:
	ClusterParam(Tcl_Interp *, const IO::ParameterManager *pPM, const IO::FileParameters * pPar);
	virtual ~ClusterParam();

	static bool   isCluster(const Kernel::ID &theMap);
	static Kernel::P_TYPE ID2pt(const Kernel::ID &);
	static Kernel::ID     pt2ID(Kernel::M_TYPE, Kernel::P_TYPE);
	static bool   isExtendedDefect(const Kernel::ID &theMap);

	static unsigned count(const std::map<std::string,bool> &array, const std::string &key);
	unsigned getDefectNumber(Kernel::M_TYPE mt, const std::string &) const;
	unsigned defectSize(Kernel::M_TYPE mt) const { return _params[mt].size(); }

	const EDType::CLType * getParams(unsigned et, const Kernel::ID &m) const;
	const EDType::CLType * getParams(Kernel::M_TYPE mt, unsigned et, unsigned hash) const;
	const EDType         * getParams(Kernel::M_TYPE mt, unsigned et) const { return _params[mt][et]; }

	std::vector<Kernel::ID> getIDs(Kernel::M_TYPE, unsigned edtype, const std::string &) const;
	bool reactionPossible(Kernel::SubDomain *, Kernel::Domain *, const Kernel::MeshElement *,
			unsigned edtype, Kernel::P_TYPE pt1, Kernel::P_TYPE pt2, float kT) const;

private:
	const IO::ParameterManager *_pPM;
	std::vector<std::vector<std::string> > _allClusters[Kernel::MAX_MATERIALS]; // mt defect_type name
	std::vector<EDType *> _params[Kernel::MAX_MATERIALS]; //storage

	void readParameters(Tcl_Interp *, const IO::ParameterManager *pPM, const IO::FileParameters * pPar);
	void readInteractions(const IO::FileParameters * pPar);
};

}

#endif /* CLUSTERPARAM_H_ */
