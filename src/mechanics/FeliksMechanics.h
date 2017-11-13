/*
 * FELIKSMechanics.h
 *
 *  Created on: May 8, 2013
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

#ifndef FELIKSMECHANICS_H_
#define FELIKSMECHANICS_H_

#include "MechInterface.h"
#include "feliks/Main/feliksinterface.h"
#include "kernel/Material.h"

class FEanalysis;

namespace Mechanics {

class FELIKSMechanics: public Mechanics::MechInterface {
public:
	FELIKSMechanics(Kernel::Domain *);
	virtual ~FELIKSMechanics();

	virtual void import() const;

private:
	std::vector<double>    				_vertices;
	std::vector<int>       				_connectivity;
	std::map<int, char>  				_bc;
	std::vector<feliks::mmoncaMaterial> 	_materials;
	mutable std::vector<feliks::elementState>   	_elstates;

	unsigned _material2idx[Kernel::MAX_MATERIALS];  //material name to idx in _materials

	unsigned idx2vert(unsigned i, unsigned j, unsigned k) const; //ijk are positions in the MMonCa grid
	unsigned _nX, _nY, _nZ;

	FEanalysis *_theAnalysis;
	feliks::mmoncaInterface *_pmi;
		//																						defect name    number of defects
	void setMechanicalState(const Kernel::MeshElement *, double &, double &, double &, std::map<std::string, unsigned> &) const;
};

} /* namespace Mechanics */
#endif /* FELIKSMECHANICS_H_ */
