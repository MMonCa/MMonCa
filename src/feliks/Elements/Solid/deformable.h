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
 * deformable.h
 *
 * geneal deformable solid
 * iro, july 2006
 */

#ifndef _deformableET_h
#define _deformableET_h

#include "Elements/eltype.h"
#include "Elements/Interpolation/interpolation.h"
#include "Elements/element.h"
#include "Math/tensor.h"
#include <iostream>

// maximum number of nodes and Gausspoints for a deformable element
#define SOLID_MAXN 10
#define SOLID_MAXG 8


class commmandLine;
class energies;
class node;



class deformableET : public eltype
{		
public:
                        deformableET(commandLine &cl);
                        deformableET();
    virtual             ~deformableET(){}
        
        
    virtual void        print(std::ostream &of=std::cout);
    virtual element*    createRandomElement() const;

};




class deformableElement : public element
	{
		
	public:
                        deformableElement(const int label, const eltype &type, 
                                          const feliks::topology::cell* c,
                                          const int nvertices, node **ndlist);
		virtual         ~deformableElement();
		
		virtual bool	check();
		virtual bool	lumpMassToNodes();
		
		bool			initialize();
		void			materialVelocityGradient(const FEshapefunbuilder& shp, blue::itensor &GradV);
        
        
        virtual bool    integrateEnergy(energies& ener, const dofset::evaluation_time& when) const=0;
        virtual bool    integrateDEnergy()=0;
		virtual bool    integrateMassMatrix(assembler& theAssembler)=0;

        
		virtual double	maxEigenvalue();
		
		
	protected:
		double			some_dimension;
	};


void			computeDeformationGradient(const std::vector<shapefun>& shp, 
                                           blue::itensor &F, 
                                           const dofset::evaluation_time when=dofset::tna);

void			pushForwardShapeFunctions(const blue::itensor &F, 
                                          const std::vector<shapefun>& shp,
                                          std::vector<blue::ivector>& b);

void			pullbackShapeFunctions(const blue::itensor &F, 
                                          const std::vector<shapefun>& shp,
                                          std::vector<blue::ivector>& b);

void            smallStrainTensor(const std::vector<shapefun>& shp,
                                  const std::vector<Unode*> nodes,
                                  const dofset::evaluation_time& when,
                                  blue::istensor &strain);

void computeIncrementalDeformationGradient(const vector<shapefun>& shp, 
                                       const blue::itensor &Fn,
                                       blue::itensor &F, 
                                       blue::itensor &g,
                                       const dofset::evaluation_time when=dofset::tna);



#endif


