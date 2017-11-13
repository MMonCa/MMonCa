/* Finite Element Method Module
 *
 * Author: ignacio.romero@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain,
 *      and       Technical University of Madrid (UPM), Madrid, Spain
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
/*
 * solid.h
 *
 * solid element infinitesima strain
 * iro, september 2002
 */

#ifndef _solid_h
#define _solid_h

#include "solidpoint.h"
#include "Model/Node/node.h"
#include "Elements/Solid/deformable.h"
#include "Math/matrix.h"
#include "Math/tensor.h"
#include "Math/vector.h"
#include "Materials/Smallstrain/alloy.h"

#include <iostream>
#include <set>

class assembler;
class commandLine;
class energies;
class quadpoint;
class resultdata;
class smallStrainMaterial;
class smallStrainMP;

namespace feliks
{
    class elementState;
};



class solid : public deformableET
{

public:
                                solid(const smallStrainMaterial& mat);
    
	virtual element*            createElement(const int label, const feliks::topology::cell* c,
                                              const int nv, node **ndlist, feliks::elementState& es) const;
    
   	virtual evalspot*           createEvalspot(const double volume, const std::vector<shapefun>& theShp, const double initspacing) const;
	virtual void                print(std::ostream &of=std::cout);
    static bool                 test();
    
private:
	bool                        mixedFormulation;
    const smallStrainMaterial*  theSmallStrainMaterial;

	friend class                element;
	friend class                solidElement;
};




class solidElement : public deformableElement
{

public:
                                solidElement(const int label, const solid &et, 
                                             const feliks::topology::cell* c,
                                             const int nvert, node **ndlist,
                                             feliks::elementState& es);
	virtual                     ~solidElement();
	virtual bool                initialize();

    virtual bool                sharedFunction(assembler& theAssembler);
    virtual bool                integrateEnergy(energies& ener, const dofset::evaluation_time& when) const;
    virtual bool                integrateDEnergy();
    virtual bool                integrateMassMatrix(assembler& theAssembler);

    virtual bool                gradient(const size_t ip, resultdata& gr, blue::ivector& gpcoor, double& dvol) const;
    
    
    // info / members
    std::set<evalspot*>         getEvalspots() const;

    
    // virtual functions to update info from time step to time step
    virtual void                commitCurrentState();
    virtual void                resetCurrentState();
	virtual void                updateCurrentState(const dofset::evaluation_time when);
    


	
private:
	const solid*                theSolidET;
	std::set<solidEvalPoint*>   theEvaluationPoints;
    feliks::elementState*         theState;
};


#endif


