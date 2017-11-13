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
	
int InitCmd::operator()()
{
	Domains::MCClient *pCli = 0;
	Kernel::Coordinates m,M;
	if(specified("mesh"))
	{
		string filename = getString("mesh");
		float scale = specified("scale")? getFloat("scale") :1000;
		IO::MeshParser dfObject(filename, scale);
		dfObject.getCellSize(m, M);
		pCli = new Domains::MCClient(&dfObject);
		Domains::global()->setClient(pCli, m, M);
		pCli->resetGetMaterial();
	}
	else
	{
		m = Coordinates(getFloat("minx"), getFloat("miny"), getFloat("minz"));
		M = Coordinates(getFloat("maxx"), getFloat("maxy"), getFloat("maxz"));
		pCli = new Domains::MCClient(_pTcl, m, M, getString("material"));
		Domains::global()->setClient(pCli, m, M);
	}

	if(specified("temp"))
		Domains::global()->setTempK(getFloat("temp") + 273.15);

	return TCL_OK;	
}


