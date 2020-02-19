/*
 * MobileParticleParam.h
 *
 *  Created on: Feb 24, 2011
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

#ifndef MOBILEPARTICLEPARAM_H_
#define MOBILEPARTICLEPARAM_H_
#include "kernel/Material.h"
#include "kernel/ParticleType.h"
#include "kernel/Event.h"
#include "kernel/Coordinates.h"
#include "io/ArrheniusAlloys.h"
#include "AlloyParam.h"

#include <map>

namespace IO { class FileParameters;  class ParameterManager; }
namespace Kernel { class MeshElement; class Domain; class StateManager; }

namespace OKMC {

class ClusterParam;

class MobileParticleParam
{
public:
	MobileParticleParam(const IO::ParameterManager *, const IO::FileParameters *);
	void init(const IO::ParameterManager *, const IO::FileParameters *, const OKMC::ClusterParam *,
				const OKMC::AlloyParam *, const Kernel::StateManager *);

	~MobileParticleParam() {}
	IO::ArrheniusAlloys _arr[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_STATES][10]; //0 mig; 1,6 br1; 2,7 br2; 3,8 FT(I); 4,9 FT(V); 5 states.
	IO::ArrheniusAlloys _orig[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_STATES][10]; //original values from the parameters
	IO::ArrheniusAlloys _form[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES]; // formation energies and prefactors	
	Kernel::Coordinates _axes[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_STATES]; //for 1D migration
	bool                _oneD[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_STATES]; //for 1D migration
	bool  				        _canInteract[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_PARTICLES];
	Kernel::Event::E_TYPE          _interact[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_PARTICLES];
	float                 _interact_radiusSq[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_PARTICLES];
	float                 _interact_radius  [Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_PARTICLES];
	std::map<unsigned, float>        _mctype[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_PARTICLES]; //which extended defects
	Kernel::P_TYPE               _int_result[Kernel::MAX_PARTICLES][Kernel::MAX_PARTICLES];

	struct MPPBreakUp
	{
		MPPBreakUp() : _isCluster(false) { _emit[0] = _emit[1] = Kernel::UNDEFINED_TYPE; }
		bool _isCluster;
		Kernel::P_TYPE _emit[2];
	} _breakUp[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_STATES];

	// parameters used for Poisson solver
	double                    _orbitalRadius[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES];
	double                    _stateEnergy[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES];     // FIXME: charge states should be used instead
	double                    _stateDegeneracy[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES]; // FIXME: charge states should be used instead
	bool                      _mapToGrid[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES];
	//
	//parameters for charge model
	bool                      _bCharged    [Kernel::MAX_MATERIALS]; //false means no charged states
	float                     _chargeLevel [Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_STATES]; //from valence band, indexed using MID_STATE!!
	int                       _state2charge[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_STATES]; //mapping.
	unsigned                  _charge2state[Kernel::MAX_MATERIALS][Kernel::MAX_PARTICLES][Kernel::MAX_STATES]; //mapping, indexed using MID_STATE for 0.

private:
	void fillResults(Kernel::M_TYPE mt, Kernel::P_TYPE pds[Kernel::MAX_PARTICLES][Kernel::MAX_PARTICLES], bool ivs[Kernel::MAX_PARTICLES][Kernel::MAX_PARTICLES]) const;
	void readReactions(const IO::ParameterManager *, const IO::FileParameters *, const OKMC::ClusterParam *);
	void fill_interact(const IO::ParameterManager *, const IO::FileParameters *);
};

}

#endif /* MOBILEPARTICLEPARAM_H_ */
