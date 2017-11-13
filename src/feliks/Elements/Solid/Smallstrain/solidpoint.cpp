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
 * solidpoint.cpp
 *
 * FELIKS,
 * Ignacio Romero, Universidad Politecnica de Madrid.
 *
 * infinitesimal strain, solid deformable element
 * It should work for 2D and 3D elements (triangles, quads, tets, bricks,...)
 *
 * i. romero, february 2002
 *            revised november 2006
 *
 */

#include "solidpoint.h"
#include "Math/tensor.h"

#include <iomanip>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include "boost/foreach.hpp"

#include "Analysis/assembler.h"
#include "Elements/Solid/deformable.h"
#include "Elements/Interpolation/interpolation.h"
#include "Geometry/quadpoint.h"
#include "General/continuum.h"
#include "Io/logger.h"
#include "Model/model.h"
#include "Model/Parts/modelpart.h"
#include "Math/matrix.h"
#include "Math/vector.h"
#include "Materials/Smallstrain/smallstrain.h"

using namespace blue;


solidEvalPoint :: solidEvalPoint(const double                   volx,
                                 const smallStrainMaterial&     mat,
                                 const std::vector<shapefun>&   sh,
                                 const feliks::elementState&      es) :
    mixedFormulation(false),
    volume(volx),
    elState(&es)
{
    theShp = sh;
    theMaterialPoint = dynamic_cast<smallStrainMP*>(mat.createMaterialPoint(&(elState->composition_)));
    
    theNodes.resize(theShp.size());
    for (size_t a=0; a<theShp.size(); a++)
        theNodes[a] = dynamic_cast<Unode*>( theShp[a].parentNode );
}




solidEvalPoint :: solidEvalPoint(const double                   volx,
                                 const smallStrainMaterial&     mat,
                                 const std::vector<shapefun>&   sh,
                                 const std::vector<shapefun>&   shBar) :
    mixedFormulation(true),
    volume(volx),
    theShpBar(shBar),
    elState(0)
{
    theShp = sh;
    theMaterialPoint = dynamic_cast<smallStrainMP*>(mat.createMaterialPoint(&(elState->composition_)));
    
    theNodes.resize(theShp.size());
    for (size_t a=0; a<theShp.size(); a++)
        theNodes[a] = dynamic_cast<Unode*>( theShp[a].parentNode );
}




solidEvalPoint :: ~solidEvalPoint()
{
    delete theMaterialPoint;
    theNodes.clear();
    theShp.clear();
    theShpBar.clear();
}




bool solidEvalPoint :: check() const
{
    bool ret = true;
        
    if ( getVolume() <= 0.0 )
    {
        logger::warnings << "\n" << "Warning: Negative jacobian in solid eval point with coordinates " << coordinates();
        ret = false;
    }
        
    return ret;
}




void solidEvalPoint :: commitCurrentState()
{
    theMaterialPoint->commitCurrentState();
}




ivector solidEvalPoint :: coordinates( const dofset::evaluation_time when) const
{
    ivector c;
    c.setZero();
    for (size_t a=0; a<theShp.size(); a++)
        c += theShp[a].value * theNodes[a]->coordinates(when);
    
    return c;
}




ivector solidEvalPoint :: getReferenceCoordinates() const
{
    ivector c;
    c.setZero();
    for (size_t a=0; a<theShp.size(); a++)
        c += theShp[a].value * theNodes[a]->getReferenceCoordinates();
    
    return c;
}




void solidEvalPoint :: computeStrain(const dofset::evaluation_time when, istensor& eps) const
{
    smallStrainTensor(this->theShp, this->theNodes, when, eps);
    double vol_strain = 0;
    // correct strain definition if mixed
    double theta=0.0;
    if (mixedFormulation)
    {
        for (size_t a=0; a<theNodes.size(); a++)
            theta += theShpBar[a].dxyz.dot( theNodes[a]->displacement(dofset::tn1) );
        
        vol_strain += (theta - eps.trace());
    }
    
    eps(0,0) -= elState->iexx_;
    eps(1,1) -= elState->ieyy_;
    eps(2,2) -= elState->iezz_;
    
    ivector brickcenter(0.0, 0.0, 0.0);
    for (size_t a=0; a<theShp.size(); a++)
        brickcenter += theNodes[a]->getReferenceCoordinates();
    brickcenter *= 1.0/theShp.size();
    
    // this is only for silicon boxes of a0 = 1.268 nm
    // best fitting with parabola of constants a and b
    const map<string, pair<double, double> > &theMap = theMaterialPoint->parentMaterial().eigenstrains_;
    double r = (getReferenceCoordinates()-brickcenter).norm();
    double L = 1.268;
   // unsigned hm = 0;
    for(map<string, pair<double, double> >::const_iterator it=theMap.begin(); it!= theMap.end(); ++it)
    {
    	map<string, unsigned>::const_iterator nDef = elState->defects_.find(it->first);
		if (nDef != elState->defects_.end()) //defect found
		{
			const double theta = -it->second.second * (r*r/L/L - it->second.first*it->second.first);
			vol_strain += nDef->second * theta;
		}
    }
    eps += vol_strain / 3.0 *istensor::identity();
}


void solidEvalPoint :: getDofNumbers(vector<int>& theDofs) const
{
    // now we fill up the idata with the degrees of freedom
    theDofs.clear();
    theDofs.reserve( 3 * theNodes.size() );
    for (size_t k=0; k<theNodes.size(); k++)
    {
        for (size_t j=0; j<3; j++)
			theDofs.push_back(theNodes[k]->getID(j));
    }
}




double solidEvalPoint :: getVolume() const
{
    return volume;
}




double solidEvalPoint :: getWaveVelocity()
{
    return 0.0;
}




bool solidEvalPoint :: gradient (resultdata& gr) const
{
    bool		computed(false);
    
    // strain at Gauss point
    smallStrainMP&  mp = *(this->theMaterialPoint);
    istensor    epsilon;
        
    istensor sigma;
    
    
    gr.setZero();
    switch ( gr.getCode() )
    {
        case RESULT_STRESS:
            mp.stress(sigma);
            gr = sigma;
            computed = true;
            break;
            
        case RESULT_STRAIN:
            computeStrain(dofset::tn1, epsilon);
            gr = epsilon;
            computed = true;
            break;
            
        case RESULT_PRESSURE:
            mp.stress(sigma);
            gr = -sigma.trace()/3.0;
            computed = true;
            break;
            
        case RESULT_PLASTIC_SLIP:
            gr = mp.plasticSlip();
            computed = true;
            break;
            
        default:
            computed = false;
    }
    
    return computed;
}




bool solidEvalPoint :: hasConstrainedDofs()
{
    bool ret = false;
    
    BOOST_FOREACH(Unode* nd, theNodes)
        ret = ret || nd->hasConstrainedDofs();
    
    return ret;
}




double solidEvalPoint :: mass() const
{
    return volume * theMaterialPoint->density();
}




bool solidEvalPoint :: genericIntegral(assembler& theAssembler)
{
    const elmtTasks& task = theAssembler.getTasks();

    // Gauss loop for numerical integration
    double rho   = theMaterialPoint->density();
    double dvol  = this->volume;
    double dmass = rho*dvol;
        
        
    ivector acc;
    if (task.dynamic_residual)
    {
        acc.setZero();
        for (size_t b=0; b< theNodes.size(); b++)
            acc += theShp[b].value * theNodes[b]->acceleration();
    }
    
    
    if (task.static_residual || task.dynamic_residual)
    {
        for (size_t a=0; a< theNodes.size(); a++)
        {
            ivector  fint; fint.setZero();
            if (task.static_residual)
            {	
                istensor sigma;
                theMaterialPoint->stress(sigma);
                sigma *= dvol;
                
                fint -= sigma*theShp[a].dxyz;            
                if (mixedFormulation)
                    fint -= sigma.trace()/3.0 * ( theShpBar[a].dxyz - theShp[a].dxyz );                
            }
            
            if (task.dynamic_residual)
                fint -= (theShp[a].value*dmass) * acc;
            
            theAssembler.assemble( fint, theNodes[a]->getUDS() );        
        }
    }    
    
    
    // Accumulate static tangent
    if (task.static_tangent)
    {
        itensor     T;
        extern  double  ftgstiff;
        double cvol = theMaterialPoint->volumetricStiffness();
        for (size_t a=0; a< theNodes.size(); a++)
            for (size_t b=0; b<theNodes.size(); b++)
            {
                if (mixedFormulation)
                {
                    theMaterialPoint->contractWithDeviatoricTangent(theShp[a].dxyz, theShp[b].dxyz, T);
                    T += cvol * itensor::dyadic(theShpBar[a].dxyz, theShpBar[b].dxyz);
                }
                else
                    theMaterialPoint->contractWithTangent(theShp[a].dxyz, theShp[b].dxyz, T);
                
                theAssembler.assemble( T*dvol*ftgstiff, theNodes[a]->getUDS(), theNodes[b]->getUDS() );
            }
    }
    
    
    // accumulate inertial tangent (consistent)
    if (task.dynamic_tangent)
    {
        extern double ftgmass;
        for(size_t a=0; a< theNodes.size(); a++)
        {
            for (size_t b=0; b< theNodes.size(); b++)
            {
                double tmp = theShp[a].value * theShp[b].value*dmass*ftgmass;
                theAssembler.assemble( tmp*itensor::identity(), theNodes[a]->getUDS(), theNodes[b]->getUDS() );
            }
        }
    }
    
    
    return true;
}




bool solidEvalPoint :: integrateDampingMatrix(assembler& a)
{
    return true;
}




bool solidEvalPoint :: integrateEnergy(energies &ener, const dofset::evaluation_time& when)
{
    // Gauss loop for numerical integration
    double rho = theMaterialPoint->density();
    double dvol  = this->volume;
    double dmass = rho*dvol;
    
     
    ivector vel, phi;
    vel.setZero();
    phi.setZero();
    for (size_t b=0; b< theNodes.size(); b++)
    {
        vel += theShp[b].value * theNodes[b]->velocity(when);
        phi += theShp[b].value * theNodes[b]->coordinates(when);
    }
    
    ener.kinetic         += .5 * vel.squaredNorm() * dmass;    
    ener.potential       += theMaterialPoint->storedEnergy()*dvol;
    ener.linearMomentum  += dmass * vel;
    ener.angularMomentum += dmass * phi.cross(vel);
     
    return true;
}




bool solidEvalPoint :: integrateDEnergy()
{
    for (size_t a=0; a< theNodes.size(); a++)
    {
        ivector  fint; fint.setZero();
        
        istensor sigma;
        theMaterialPoint->stress(sigma);
        
        fint = sigma*theShp[a].dxyz;            
        if (mixedFormulation)
            fint += sigma.trace()/3.0 * ( theShpBar[a].dxyz - theShp[a].dxyz );
        
        theNodes[a]->assembleForce(fint*this->volume);
    }
    return true;
}




bool solidEvalPoint :: integrateNodalMasses()
{	
    double  rho = theMaterialPoint->density();
	if (rho <= 0.0) return false;

    double mass = rho*volume;
    
    for(size_t a=0; a<theNodes.size(); a++)
        theNodes[a]->incrementMass( theShp[a].value * mass );
    
	return true;
}




bool solidEvalPoint :: integrateMassMatrix(assembler& theAssembler)
{	
    double  rho = theMaterialPoint->density();
	if (rho <= 0.0) return false;
    
    double mass = rho*volume;
    
    for(size_t a=0; a<theNodes.size(); a++)
    {
        for (size_t b=0; b<theNodes.size(); b++)
        {
            double mab = theShp[a].value * theShp[b].value * mass;
            
            theAssembler.assemble(mab, theNodes[a]->getUDS(), theNodes[b]->getUDS() );
        }
    }
    
	return true;
}




void solidEvalPoint :: resetCurrentState()
{
    theMaterialPoint->resetCurrentState();
}




void solidEvalPoint :: updateCurrentState(const dofset::evaluation_time when)
{
    istensor    epsilon;
    computeStrain(when, epsilon);
    
    extern double global_tna, global_tn1;
    double theTime = (when == dofset::tn1) ? global_tn1 : global_tna;
    theMaterialPoint->updateCurrentState(theTime, epsilon, elState->temperature_);
}




