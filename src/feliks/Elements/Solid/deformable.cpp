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
 * deformable.cpp
 *
 * FELIKS,
 * Ignacio Romero, Universidad Politecnica de Madrid.
 *
 * general purpose deformable solid
 * i. romero, july 2006
 *
 */

#include "Elements/Solid/deformable.h"
#include "Elements/Interpolation/interpolation.h"

#include <cstring>
#include <cmath>

#include "Analysis/assembler.h"
#include "Elements/eltype.h"
#include "Geometry/quadpoint.h"
#include "Model/model.h"
#include "Materials/material.h"
#include "Io/message.h"
#include "Math/matrix.h"
#include "Math/vector.h"
#include "Io/logger.h"

using namespace blue;


deformableET :: deformableET(commandLine &cl) :
    eltype(cl)
{
	// constants
	_geometry       = VOLUME;    
}



deformableET :: deformableET() :
    eltype()
{
    _geometry = VOLUME;
}




element* deformableET :: createRandomElement() const
{
    return 0;
}




void deformableET :: print(std::ostream &of)
{
	eltype::print(of);
}



deformableElement :: deformableElement(const int label, const eltype &type, 
                                       const feliks::topology::cell* c,
                                       const int nvertices, node **ndlist) :
    element(label, type, c, nvertices, ndlist),
    some_dimension(0.0)
{}


deformableElement :: ~deformableElement()
{
}


bool deformableElement :: check()
{
	if (getNNodes() > SOLID_MAXN)
	{
		logger::warnings << endl << "Number of adminissible nodes reached in element: " << getLabel() << flush;
		return false;
	}
    
    return true;
}



bool deformableElement :: initialize()
{	
	some_dimension = characteristicDim();
	return 	element::initialize();
}



bool deformableElement :: lumpMassToNodes()
{	
	quadpoint   gauss[SOLID_MAXG];
	FEshapefunbuilder   shp(*this);	
	
    // Calculate Gauss points and weights for quadrature
	size_t     ngp;
	double  jacobian;
    fullQuadraturePoints(getNNodes(), VOLUME, &ngp, gauss);
	
    // Get material data
    const material *mat = &getMaterial();
	
	if ( mat != NULL)
	{
		double rho    = mat->density();
		
		// Gauss loop for numerical integration
		for (int ip=0; ip<ngp; ip++)
		{
			shp.evaluate(gauss[ip], jacobian);
            
            if (jacobian < 0.0) {
                jacobian = -jacobian;        //MMM Not clear whether this should be done here.
            }
			double dmass = rho * jacobian * gauss[ip].getWeight();
			
			for(int a=0; a<getNNodes(); a++)
				getNode(a).incrementMass( shp[a].value*dmass );
		}
	}
	lumpedMass=true;
	return true;
}



/* this is the function that computes the different mass matrices. Notice that it is
 * common for all deformable elements
 
bool deformableElement :: integrateMassMatrix(assembler& theAssembler)
{
	quadpoint   gauss[SOLID_MAXG];
	FEshapefunbuilder   shp(*this);	
    
	
    // Calculate Gauss points and weights for quadrature
	int ngp;
	double jacobian;
    fullQuadraturePoints(getNNodes(), VOLUME, &ngp, gauss);
	
    // Get material data
    const material& mat = getMaterial();
	istensor Mab;
	
    double rho = mat.density();
    
    
    // Gauss loop for numerical integration
    for (int ip=0; ip<ngp; ip++)
    {
        shp.evaluate(gauss[ip], jacobian);
        double dvol  = jacobian * gauss[ip].getWeight();
        double dmass = rho*dvol;
        
        for(int a=0; a<getNNodes(); a++)
        {
            for (int b=0; b<getNNodes(); b++)
            {
                double tmp = shp[a].value*shp[b].value*dmass;
                Mab = tmp * istensor::identity();
                
                theAssembler.assemble(Mab, dynamic_cast<Unode&>(getNode(a)).getUDS(), dynamic_cast<Unode&>(getNode(b)).getUDS()  );
            }
        }
    }
	return true;
}
*/


void deformableElement :: materialVelocityGradient(const FEshapefunbuilder& shp, itensor &GradV)
{
	GradV.setZero();
	for (int a=0; a<shp.size(); a++)
	{		
		// GradV += Va otimes Grad Na
		GradV.addDyadic( getNode(a).velocity(), shp[a].dxyz );
	}	
}




double deformableElement :: maxEigenvalue()
{
	return 0.0;
}




void computeDeformationGradient(const vector<shapefun>& shp, itensor &F, const dofset::evaluation_time when)
{
	F = itensor::identity();
	
	for (size_t a=0; a< shp.size(); a++)
		F.addDyadic( shp[a].parentNode->getUDS().displacement(when), shp[a].dxyz );
}




void computeIncrementalDeformationGradient(const vector<shapefun>& shp, const itensor &Fn, itensor &F, itensor &g, const dofset::evaluation_time when)
{
	g = itensor::identity();
    
    const size_t nn = shp.size();
	for (size_t a=0; a<nn; a++)
    {
        node* nd = shp[a].parentNode;
        
        ivector deltau = nd->getUDS().displacement(when);
        deltau -= nd->getUDS().displacement(dofset::tn);
        
        g.addDyadic( deltau, shp[a].dxyz );
    }
    F = g * Fn;
}




void  pushForwardShapeFunctions(const itensor &F, const vector<shapefun>& shp, vector<ivector>& b)
{
	itensor iF = F.inverse();
	itensor iFt = iF.transpose();
    b.clear();
	
    for (int a=0; a< shp.size(); a++)
		b.push_back( iFt * shp[a].dxyz ); 
}




void  pullbackShapeFunctions(const itensor &F, const vector<shapefun>& shp, vector<ivector>& b)
{
	itensor Ft = F.transpose();
    b.resize(shp.size());
	
    for (int a=0; a< shp.size(); a++)
		b[a] = Ft * shp[a].dxyz;
}




/*  compute infinitesimal strain for 3d solid.
 */
void smallStrainTensor(const vector<shapefun>& shp, 
                       const vector<Unode*> nodes, 
                       const dofset::evaluation_time& when,
                       istensor &strain)
{
	itensor beta;
    beta.setZero();

    for (size_t a=0; a<shp.size(); a++)
		beta.addDyadic( nodes[a]->getUDS().displacement(when), shp[a].dxyz);
    
    strain = istensor::symmetricPartOf(beta);
}




