/*
 * ProfileCmd.cpp
 *
 *  Created on: Mar 10, 2011
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

#include "ProfileCmd.h"
#include "kernel/ParticleType.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include "kernel/RateManager.h"
#include "domains/MCClient.h"
#include "io/ParameterManager.h"
#include "okmc/MobileParticle.h"
#include "okmc/Cluster.h"
#include "okmc/AlloyParam.h"
#include "lkmc/LKMCModel.h"
#include <sstream>
#include <cstdlib>

using std::stringstream;
using std::string;
using Kernel::Coordinates;

namespace IO {

ProfileCmd::ProfileCmd(Tcl_Interp *p, int argc, const char *argv[]) : Command(p, argc, argv)
{
}

int ProfileCmd::operator()()
{
	Domains::global()->client()->beginInsert();
	if(specified("proc"))
	{
		string proc = getString("proc");
		string name = getString("name");
		string stat = specified("state")? getString("state") : "";
		bool bReact = !specified("do.not.react");
		for(Domains::MeshElementIterator it = Domains::global()->beginMEI(); it != Domains::global()->endMEI(); ++it)
		{
			Coordinates m, M;
			it->getCorners(m, M);
			Coordinates c;
			(*it)->getDomain()->_pMesh->getCenter(c,(*it)->getIndex());
			//run the script...
			stringstream cmd;
			cmd << proc << " " << c._x << " " << c._y << " " << c._z;
			if(Tcl_EvalEx(_pTcl,  cmd.str().c_str(), -1, 0) != TCL_OK)
				ERRORMSG("The execution of the profile script named '" << proc << "' failed.");
			//ret is in cm^-3 for historical reasons
			string ret  = Tcl_GetStringResult(_pTcl);
			double conc = std::atof(ret.c_str());

			Kernel::M_TYPE mt = it->getMaterial();
			if(conc * 0.99 > Domains::global()->PM()->getMaterial(mt)._densityAlloyCm3 &&
					Domains::global()->PM()->getMaterial(mt)._alloyName == name)
				ERRORMSG("Alloy profile must be equal or lower than alloy concentration [" <<
						Domains::global()->PM()->getMaterial(mt)._densityAlloyCm3 << "]");

			// I have a conc in a box m, M with volume
			float number_part = it->getVolume() * 1e-21 * conc;
			Kernel::Event::E_TYPE et = Domains::global()->PM()->getDefectType(name, mt);

			while(Domains::global()->client()->rand() < number_part--)
			{
				Coordinates mpcoord(
						m._x + (M._x - m._x)*Domains::global()->client()->rand(),
						m._y + (M._y - m._y)*Domains::global()->client()->rand(),
						m._z + (M._z - m._z)*Domains::global()->client()->rand());

				switch(et)
				{
				case Kernel::Event::MOBILEPARTICLE:
					Domains::global()->client()->createMP(name, stat, mpcoord, bReact);
				break;
				case Kernel::Event::CLUSTER:
				{
					string defect = getString("defect");
					Domains::global()->client()->createMC(name, defect, mpcoord);
				}
				break;
				default:
					WARNINGMSG("Sorry, defect '" << name << "' not implemented yet in the material " << Domains::global()->PM()->getMaterialName(mt));

					if(Domains::global()->PM()->getMaterial(mt)._bAmorphous)
						if(name == "I" || name == "V" )
							(*it)->getDomain()->_pMesh->getElement(it->getIndex())->incrAmorphParts();

					break;
				}
			}
			//after having "inserted" an alloy, the LatticeAtoms should be checked!
			if(Domains::global()->PM()->getMaterial(mt)._alloyName == name)
			{
				Kernel::P_TYPE alloy = Domains::global()->PM()->getMaterial(mt)._alloy;
				Kernel::P_TYPE matPt =  Domains::global()->PM()->getMaterial(mt)._pt[0];
				float prob = (*it)->getBAtoms() / float((*it)->getAtoms());
				LKMC::LatticeSite *pLS = (*it)->getFirstLS();
				while(pLS)
				{
					if(pLS->getPType() == alloy || pLS->getPType() == matPt)
					{
						Kernel::P_TYPE myPt =  (Domains::global()->client()->rand() < prob? alloy : matPt);
						pLS->setPType(myPt);
					}
					//out of the if to update also neighbors indirectly...
					LKMC::LatticeAtom *pLA = dynamic_cast<LKMC::LatticeAtom *>(pLS);
					if(pLA) //no need to update neighbors, because there is a for loop updating all the cells.
						pLA->getDomain()->_pRM->update(pLA, pLA->getElement());
					pLS = pLS->getNext();
				}
			}
		}
	}
	Domains::global()->client()->endInsert();
	return TCL_OK;
}


}
