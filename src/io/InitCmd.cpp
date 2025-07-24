/*
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

#include "InitCmd.h"
#include <tcl.h>
#include "kernel/Coordinates.h"
#include "kernel/Domain.h"
#include "domains/MCClient.h"
using namespace IO;
using std::string;
using Kernel::Coordinates;

InitCmd::InitCmd(Tcl_Interp *p, int argc, const char *argv[]) : Command(p, argc, argv)
{
}

void checkOrderingUniqueness(std::vector<float> const& aVector) {
    bool was = false;
    float previous;
    for(auto const value : aVector) {
        if(was) {
            if(previous >= value) {
                ERRORMSG("InitCmd: lines* must be strictly ascending array.");
            }
        }
        else {
            was = true;
        }
        previous = value;
    }
}
	
int InitCmd::operator()()
{
	Domains::MCClient *pCli = 0;
        if(specified("mesh")) {
            std::unique_ptr<MeshMaterialReader> meshMaterial(new MeshMaterialReader(getString("mesh")));
            auto const linesX = meshMaterial->getLinesX();
            auto const linesY = meshMaterial->getLinesY();
            auto const linesZ = meshMaterial->getLinesZ();
            Kernel::Coordinates const m = Coordinates(linesX.front(), linesY.front(), linesZ.front());
            Kernel::Coordinates const M = Coordinates(linesX.back(), linesY.back(), linesZ.back());
            pCli = new Domains::MCClient(_pTcl, m, M, std::move(meshMaterial));
            Domains::global()->setClient(pCli, m, M, &linesX, &linesY, &linesZ);
        }
        else if(specified("linesx") && specified("linesy") && specified("linesz")) {
            auto const linesX = getFloats("linesx");
            auto const linesY = getFloats("linesy");
            auto const linesZ = getFloats("linesz");
            checkOrderingUniqueness(linesX);
            checkOrderingUniqueness(linesY);
            checkOrderingUniqueness(linesZ);
            Kernel::Coordinates const m = Coordinates(linesX.front(), linesY.front(), linesZ.front());
            Kernel::Coordinates const M = Coordinates(linesX.back(), linesY.back(), linesZ.back());
            pCli = new Domains::MCClient(_pTcl, m, M, getString("material"));
            Domains::global()->setClient(pCli, m, M, &linesX, &linesY, &linesZ);
        }
        else {
            Kernel::Coordinates const m = Coordinates(getFloat("minx"), getFloat("miny"), getFloat("minz"));
            Kernel::Coordinates const M = Coordinates(getFloat("maxx"), getFloat("maxy"), getFloat("maxz"));
            pCli = new Domains::MCClient(_pTcl, m, M, getString("material"));
            Domains::global()->setClient(pCli, m, M);
        }

	if(specified("temp"))
		Domains::global()->setTempK(getFloat("temp") + 273.15);
	if(specified("totalseconds"))
		Domains::global()->setTotalAllowedSeconds(getFloat("totalseconds"));

	return TCL_OK;	
}


