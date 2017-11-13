/*
 * InterfaceParam.h
 *
 *  Created on: Feb 28, 2011
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

#ifndef INTERFACEPARAM_H_
#define INTERFACEPARAM_H_

#include "io/ArrheniusAlloys.h"
#include "kernel/Material.h"
#include "kernel/ParticleType.h"
#include <vector>
#include <map>

namespace IO { class ParameterManager; class FileParameters; }
namespace Kernel { class MeshParam; }

namespace OKMC {

class MobileParticleParam;
class ClusterParam;

class InterfaceParam
{
public:
	InterfaceParam(const IO::ParameterManager *_pPM, const IO::FileParameters * pPar,
			OKMC::MobileParticleParam * pMPPar, const OKMC::ClusterParam *, const Kernel::MeshParam *);
	~InterfaceParam();

	// An InterfaceParamSet stores the parameter for a particular side of the interface.
	struct InterfaceParamSet
	{
		InterfaceParamSet(const OKMC::ClusterParam *pEDPar, Kernel::M_TYPE mt);

        IO::ArrheniusAlloys      _arrEmitMP[Kernel::MAX_PARTICLES];
		float                    _barrierMP[Kernel::MAX_PARTICLES];
		IO::Arrhenius            _formation[Kernel::MAX_IMPURITIES];
		IO::Arrhenius            _migration[Kernel::MAX_IMPURITIES];
		float					 _lambda;  //migration distance
		float					 _superSat [Kernel::MAX_PARTICLES];

		float                    _trapMPProb[Kernel::MAX_PARTICLES]; //probability for mobile particles to be trapped by the Interface
		std::vector<float>       _trapMCProb;

		float                    _desorptionMPProbL[Kernel::MAX_PARTICLES]; //probability for mobile particles to be annihilated
		float                    _desorptionMPProbH[Kernel::MAX_PARTICLES]; //range 2 for desorption
		float                    _desorptionThreshold[Kernel::MAX_PARTICLES]; //treshold between probs in cm^-2
		std::vector<float>       _desorptionMCProb;

		bool                     _interactWithMP[Kernel::MAX_PARTICLES];
		std::vector<bool>        _interactWithMC;
	};

	//Stores pointers to ALL the interface sides. First is FROM, second is TO.
	//0 means such side is not defined and should not exist in the simulation.
	InterfaceParamSet *_params[Kernel::MAX_MATERIALS][Kernel::MAX_MATERIALS]; //from -> to

private:
	unsigned _nMat;

	std::string getInterfaceName(const IO::ParameterManager *pPM, Kernel::M_TYPE mt1, Kernel::M_TYPE mt2) const;
	std::string getTo  (const IO::ParameterManager *pPM, const IO::FileParameters * pPar, const std::string &base, Kernel::M_TYPE mtTo) const;

	void readParameters(const IO::ParameterManager *pPM, const IO::FileParameters *pPar,
			OKMC::MobileParticleParam * pMPPar, const OKMC::ClusterParam *pMCPar, const Kernel::MeshParam *);
	void readReactions (const IO::ParameterManager *pPM, const IO::FileParameters *pPar,
			const OKMC::ClusterParam *pMCPar);

	InterfaceParam(const InterfaceParam &); //forbidden
	InterfaceParam & operator=(const InterfaceParam &);
};

}

#endif /* INTERFACEPARAM_H_ */
