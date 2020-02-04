/*
 * GetMaterial.cpp
 *
 *  Created on: Jun 20, 2012
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


#include "GetMaterial.h"
#include "Diagnostic.h"
#include "ParameterManager.h"
#include "domains/Splitter.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include <sstream>

using std::stringstream;
using std::string;

namespace IO {

GetMaterial::GetMaterial(std::istream &is) : _pTcl(0), _procName("")
{
	int nDomains;
	is >> nDomains;
	for(int dom=0; dom<nDomains; ++dom)
	{
		int nElements;
		is >> nElements;
		_materials.push_back(std::vector<Kernel::M_TYPE>());
		_materials.back().reserve(nElements);
		for(int i=0; i < nElements; i++)
		{
			Kernel::M_TYPE mt;
			is >> mt;
			_materials[dom].push_back(mt);
		}
	}

}

void GetMaterial::restart(std::ostream &os)
{
	int nDomains = Domains::global()->getSplitter()->getDomains();
	os << nDomains << std::endl;
	for(int dom=0; dom<nDomains; ++dom)
	{
		int nElements = Domains::global()->getDomain(dom)->_pMesh->size();
		os << nElements << std::endl;
		for(int i=0; i<nElements; ++i)
			os << Domains::global()->getDomain(dom)->_pMesh->getElement(i)->getMaterial();
	}
}

Kernel::M_TYPE GetMaterial::operator()(const Kernel::Coordinates &c, unsigned idx) const
{
	if(_materials.size())
	{
		int domain = Domains::global()->getSplitter()->getDomain(c);
		return _materials[domain][idx];
	}

	//run a tcl procedure and extract the info from it.
	stringstream cmd;
	cmd << _procName << " " << c._x << " " << c._y << " " << c._z;
	if(Tcl_EvalEx(_pTcl,  cmd.str().c_str(), -1, 0) != TCL_OK)
		ERRORMSG("The execution of the material script named '" << _procName << "' failed.");
	string mat = Tcl_GetStringResult(_pTcl);
	Kernel::M_TYPE mt = Domains::global()->PM()->getMaterialNumber(mat);
	if(mt == Kernel::MAX_MATERIALS)
		ERRORMSG("Material " << mat << " not recognized!");
	return mt;
}

} /* namespace IO */
