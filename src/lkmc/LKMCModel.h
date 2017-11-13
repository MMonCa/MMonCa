/*
 * LKMCModel.h
 *
 *  Created on: May 11, 2011
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

#ifndef LKMCMODEL_H_
#define LKMCMODEL_H_

#include <vector>
#include "LKMCMode.h"
#include "kernel/Material.h"

namespace Kernel { class Domain; class SubDomain; class MeshElement; class SubDomain; }

namespace LKMC {

class LKMCModel {
public:
	LKMCModel(bool bFromStream, Kernel::Domain *);
	virtual ~LKMCModel() {}
	unsigned putLKMCAtoms(Kernel::SubDomain *, Kernel::MeshElement *, Kernel::M_TYPE basicMat, LKMCMode) const;
	void cleanLKMCAtoms(Kernel::MeshElement *pEle, LKMCMode) const;
	void localSPER(Kernel::SubDomain *, const std::vector<Kernel::MeshElement *> &elems) const;
private:
	Kernel::Domain *_pDomain;
	void epitaxy() const;
	void sper() const;

};

}

#endif /* LKMCMODEL_H_ */
