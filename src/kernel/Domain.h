/* Author: ignacio.martin@imdea.org
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

#ifndef KMCDOMAIN_H
#define KMCDOMAIN_H

#include "Coordinates.h"
#include "Material.h"
#include "RNG.h"
#include "UpdateManager.h"

#include <string>
#include <vector>

namespace LKMC
{
	class Lattice;
	class LatticeParam;
	class LKMCModel;
	struct EpiGasParam;
}

namespace OKMC
{
	class MobileParticleParam;
	class AlloyParam;
	class InterfaceParam;
	class ClusterParam;
	class ClusterReactionParam;
	class OKMCModel;
}

namespace Mechanics { class MechInterface; }

namespace IO
{
	class FileParameters;
	class Global;
}

namespace Electrostatics {
    class Poisson;
}

namespace Domains {
	class MCClient;
}

struct Tcl_Interp;

namespace Kernel
{
class Mesh;
class RateManager;
class MeshParam;
class StateManager;

class Domain
{
public:
    Domain(unsigned num, Tcl_Interp *, const Coordinates &m, const Coordinates &M,
           std::vector<float> const * const aLinesX = nullptr, std::vector<float> const * const aLinesY = nullptr, std::vector<float> const * const aLinesZ = nullptr);
    void init(bool bFromStream); //is the simulation recreated from a stream?
    ~Domain();
    
    const LKMC::LatticeParam         *_pLaPar[MAX_MATERIALS];
    	  LKMC::EpiGasParam          *_pEGPar[MAX_MATERIALS];

    const MeshParam                  *_pMePar;
          OKMC::AlloyParam           *_pAlPar;
    const OKMC::MobileParticleParam  *_pMPPar;
    const OKMC::InterfaceParam       *_pIFPar;
    const OKMC::ClusterParam    	 *_pClPar;
    const OKMC::ClusterReactionParam *_pClRePar;

	Mesh           *_pMesh;
	RateManager    *_pRM;

        Electrostatics::Poisson        *_pPoisson;
	const LKMC::Lattice  *_pLat[MAX_MATERIALS];
	
	const LKMC::LKMCModel *_pLKMC;
	const OKMC::OKMCModel	*_pOKMC;

	StateManager		         *_pSM;
	const Mechanics::MechInterface * _pMech;

	RNG _rng_dom;  //for MC
	UpdateManager _um_mechani;
	UpdateManager _um_poisson;
	const unsigned _domain_number;
};

}

#endif
