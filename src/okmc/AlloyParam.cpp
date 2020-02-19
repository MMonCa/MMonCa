/*
 * AlloyParam.cpp
 *
 *  Created on: Mar 11, 2013
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

#include "okmc/AlloyParam.h"
#include "Defect.h"
#include "io/FileParameters.h"
#include "io/ParameterManager.h"

using IO::ArrheniusAlloys;
using Kernel::M_TYPE;
using Kernel::P_TYPE;
using std::map;
using std::string;
using std::vector;

namespace OKMC {

AlloyParam::AlloyParam(const IO::ParameterManager *pPM, const IO::FileParameters * pPar)
{
	for(M_TYPE mt=0; mt < pPM->getNMaterials(); ++mt)
        {       
                _isAlloy[mt] = pPM->getMaterial(mt)._alloyName != "none" ? true : false;
		if(_isAlloy[mt])
		{
			if(pPM->getMaterial(mt)._binary)
				ERRORMSG("Alloying not implemented yet for binary materials");
			std::stringstream ss_selfDif;
			ss_selfDif << pPM->getMaterialName(mt) << "/Models/self.diffusion";
			if(!pPar->getBool(ss_selfDif.str()))
				WARNINGMSG("Alloy material is given but selfdiffusion is set to false");
			std::stringstream ss_theta;
			ss_theta << pPM->getMaterialName(mt) << "/Models/theta";
			_theta[mt] = pPar->getFloat(ss_theta.str());
			_mixingEnthalpy[mt] = pPar->getPolynomial(pPM->getMaterialName(mt) + "/Models/mixing.enthalpy");
			float L = pPar->getFloat("MC/Mesh/spacing.x");
			float r = pPar->getFloat (pPM->getMaterialName(mt) + "/Lattice/" + "parameter") / pPar->getFloat("MC/Mesh/spacing.x");
			float w2th = pPar->getFloat (pPM->getMaterialName(mt) + "/Models/alloy.second.neighbor.contribution");
			float f = pPar->getFloat (pPM->getMaterialName(mt) + "/Models/smoothing.factor");
			if(L != pPar->getFloat("MC/Mesh/spacing.y") || L != pPar->getFloat("MC/Mesh/spacing.z"))
				WARNINGMSG("Alloy concentration smoothing algorithm requires cubic Mesh Elements.");
			if(SVAL != 3)
				ERRORMSG("Alloy concentration smoothing algorithm expects 26 neighbors");
			_smoothing[mt][0] = 3 * r / (8 + 6 * w2th) * (r * r - 4 * r + 4 + 2 * w2th);
			_smoothing[mt][1] = 3 * r * r / (8 + 6 * w2th) * (2 - r);
			_smoothing[mt][2] = r * r * r / (8 + 6 * w2th);

			MEDMSG("Concentration smoothing algorithm parameters for " << pPM->getMaterialName(mt));
			MEDMSG("Smoothing factor: " << f);
			MEDMSG("Second neighbor contribution: " << w2th);
			MEDMSG("Relative weight for the face boxes " << _smoothing[mt][0]);
			MEDMSG("Relative weight for the edge boxes " << _smoothing[mt][1]);
			MEDMSG("Relative weight for the corner boxes " << _smoothing[mt][2]);
			MEDMSG("");

			_smoothing[mt][0] *= f / 6.;   // 6 side boxes
			_smoothing[mt][1] *= f / 12.;  // 12 edge boxes
			_smoothing[mt][2] *= f / 8.;   // 8 corner boxes

			for(P_TYPE iorv=0; iorv < 2; ++iorv)
			{
				P_TYPE pt = pPM->iorv2pt(mt, iorv, Kernel::POS_0);
				std::stringstream ss_fact, ss_alpha;
				string familyName = pPM->getFamilyName(pPM->getFamily(pt));

				ss_fact  << pPM->getMaterialName(mt) << "/" << familyName << "/correlation.factor";
				ss_alpha  << pPM->getMaterialName(mt) << "/" << familyName << "/alpha";				
				_corrFactor[mt][iorv] = pPar->getFloat(ss_fact.str());
				_alpha[mt][iorv] = pPar->getArrheniusAlloys(ss_alpha.str());				
			}
		}
        }
}

double AlloyParam::getDerMixingEnergy(M_TYPE mt, double x, float T)
{	    
	return _mixingEnthalpy[mt].getIntFirstDerValue(x) * (1 - T / _theta[mt]);
}

double AlloyParam::getMixingEnergy(M_TYPE mt, double x, float T)
{
	
	return _mixingEnthalpy[mt].getIntValue(x) * (1 - T / _theta[mt]);
}

}


