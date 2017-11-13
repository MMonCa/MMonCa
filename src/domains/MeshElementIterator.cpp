/*
 * MeshElementIterator.cpp
 *
 *  Created on: Jun 20, 2012
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

#include "MeshElementIterator.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include "Global.h"

namespace Domains {


MeshElementIterator::MeshElementIterator(const MeshElementIterator &m)
{
	_pME = m._pME;
	_nDomain = m._nDomain;
	_nElement = m._nElement;
}

MeshElementIterator & MeshElementIterator::operator=(const MeshElementIterator &m)
{
	_pME = m._pME;
	_nDomain = m._nDomain;
	_nElement = m._nElement;
	return *this;
}

MeshElementIterator & MeshElementIterator::operator++()
{
	Kernel::Domain * pDomain = global()->getDomain(_nDomain);
	if(++_nElement == pDomain->_pMesh->size())  //next domain
	{
		_nElement = 0;
		if(++_nDomain < global()->getDomains())
			pDomain = global()->getDomain(_nDomain);
		else
		{
			_pME = 0;
			return *this; //with end() value, defined as _nDomain pointing out
		}
	}
	_pME = pDomain->_pMesh->getElement(_nElement);
	return *this;
}

} /* namespace Domains */
