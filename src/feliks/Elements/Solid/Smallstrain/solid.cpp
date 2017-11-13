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
 * solid.cpp
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

#include "solid.h"

#include <iomanip>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include "boost/foreach.hpp"

#include "Main/feliksinterface.h"
#include "Analysis/assembler.h"
#include "Model/model.h"
#include "Model/Parts/modelpart.h"
#include "Geometry/quadpoint.h"
#include "Elements/Interpolation/interpolation.h"
#include "Math/feliksmath.h"
#include "Io/io.h"
#include "Materials/Smallstrain/smallstrainlib.h"

using namespace blue;




solid :: solid(const smallStrainMaterial& am)
:
    deformableET(),
    mixedFormulation(false),
    theSmallStrainMaterial(&am)
{
    _name     = "SOLID";
	_nodetype = _Unode;
	_geometry = VOLUME;    
}




element* solid :: createElement(const int label, const feliks::topology::cell* c,
                                const int nv, node **ndlist, feliks::elementState& es) const
{
	return new solidElement(label, *this, c, nv, ndlist, es);
}	




evalspot* solid :: createEvalspot(const double volume, const std::vector<shapefun>& theShp, const double initspacing) const
{
    return 0;//new solidEvalPoint(volume, *theSmallStrainMaterial, theShp);
}




void solid :: print(std::ostream &of)
{
	deformableET::print(of);
    of << "\n Infinitesimal strain deformable element.";
    if (mixedFormulation)    of << "\n Mixed -- mean volume.";
    else                     of << "\n Displacement based";
}




bool solid :: test()
{
    return false;
}




solidElement :: solidElement(const int labelx, const solid &et, const feliks::topology::cell* c,
                             const int nvertices, node **ndlist, feliks::elementState& es) :
    deformableElement(labelx, et, c, nvertices, ndlist),
    theSolidET(&et),
    theState(&es)
{
	
}




solidElement :: ~solidElement()
{
    BOOST_FOREACH(solidEvalPoint* ev, theEvaluationPoints) 
        delete ev;
    theEvaluationPoints.clear();
}





void solidElement :: commitCurrentState()
{
    BOOST_FOREACH(evalspot* ev, theEvaluationPoints) ev->commitCurrentState();
}




std::set<evalspot*> solidElement :: getEvalspots() const
{
    std::set<evalspot*> theEV;
    
    BOOST_FOREACH(solidEvalPoint* ev, theEvaluationPoints)
    {
        theEV.insert(ev);
    }
    
    return theEV;
}




bool solidElement :: gradient(const size_t ip, resultdata& gr, ivector& gpcoor, double& dvol) const
{
    std::set<solidEvalPoint*>::const_iterator iter = theEvaluationPoints.begin();
    for (size_t a=0; a<ip; a++) ++iter;
    
    (*iter)->gradient(gr);
    gpcoor = (*iter)->coordinates();
    dvol   = (*iter)->getVolume();
    
    return true;
}




bool solidElement :: initialize()
{
	deformableElement::initialize();
    	
    // location of quadrature points
    size_t      ngp=0;
    quadpoint   gauss[SOLID_MAXG];
    
    switch (getNNodes()) 
    {
        case 4:     ngp = 1; break;
        case 8:     ngp = 8; break;
        case 10:    ngp = 5; break;
        default:    
            logger::warnings << "\n Element with bad number of nodes " << getLabel();
            break;
    }
    quadraturePoints(getNNodes(), VOLUME, ngp, gauss);
    
    
    // averaged shape functions
    std::vector<shapefun>  theAveragedShapeFunctions;
    if (theSolidET->mixedFormulation)
    {
        theAveragedShapeFunctions.resize( getNNodes() );
        for (size_t a=0; a<getNNodes(); a++)
        {
            theAveragedShapeFunctions[a].value = 0.0;
            theAveragedShapeFunctions[a].dxyz.setZero();
        }
        
        double vol = 0.0;
        double dvol;    
        for (size_t ip=0; ip<ngp; ip++)
        {
            FEshapefunbuilder theShapeFunctions(*this);
            double jacobian;
            theShapeFunctions.evaluate(gauss[ip], jacobian);
            dvol = jacobian*gauss[ip].getWeight();
            vol += dvol;
            
            for (size_t a=0; a< getNNodes(); a++)
            {
                theAveragedShapeFunctions[a].value += theShapeFunctions[a].value * dvol;
                theAveragedShapeFunctions[a].dxyz  += theShapeFunctions[a].dxyz  * dvol;
            }
        }	
        
        
        for (size_t a=0; a< getNNodes(); a++)
        {
            theAveragedShapeFunctions[a].value *= 1.0/vol;
            theAveragedShapeFunctions[a].dxyz  *= 1.0/vol;
        }
    }
    
    // regular shape functions and creation of evalspots
    for (size_t ip=0; ip<ngp; ip++)
    {
        FEshapefunbuilder   theShapeFunctions(*this);
        double jacobian;
        theShapeFunctions.evaluate(gauss[ip], jacobian);
        double dvol = jacobian*gauss[ip].getWeight();
        
        // revert inverted elements
        if (dvol < 0.0)
        {
            logger::warnings << "\n Inverting element " << getLabel();
            dvol = -dvol;
        }
        
        std::vector<shapefun> theShp;
        for (size_t a=0; a<getNNodes(); a++)
            theShp.push_back( theShapeFunctions[a] );
        
        solidEvalPoint* ep;
        if (theSolidET->mixedFormulation)
             ep = new solidEvalPoint(dvol, 
                                    *(theSolidET->theSmallStrainMaterial),
                                    theShp, 
                                    theAveragedShapeFunctions);
        else
            ep = new solidEvalPoint(dvol, 
                                    *(theSolidET->theSmallStrainMaterial),
                                    theShp,
                                    *theState);
            
        theEvaluationPoints.insert(ep);
    }    
	
	return true;
}




bool solidElement :: integrateEnergy(energies& ener, const dofset::evaluation_time& when) const
{
    BOOST_FOREACH(evalspot* ev, theEvaluationPoints) 
        ev->integrateEnergy(ener, when);
    return true;    
}




bool solidElement :: integrateDEnergy()
{
    BOOST_FOREACH(evalspot* ev, theEvaluationPoints) 
        ev->integrateDEnergy();
    return true;
}




bool solidElement :: integrateMassMatrix(assembler& theAssembler)
{
    BOOST_FOREACH(evalspot* ev, theEvaluationPoints)
        ev->integrateMassMatrix(theAssembler);
    
    return true;
}




void solidElement :: resetCurrentState()
{
    BOOST_FOREACH(evalspot* ev, theEvaluationPoints) ev->resetCurrentState();
}




bool solidElement :: sharedFunction(assembler& theAssembler)
{
    BOOST_FOREACH(evalspot* ev, theEvaluationPoints)
        ev->genericIntegral(theAssembler);
    return true;
}




void solidElement :: updateCurrentState(const dofset::evaluation_time when)
{
	istensor strain, stress;
	strain.setZero();
    stress.setZero();

    BOOST_FOREACH(solidEvalPoint* ev, theEvaluationPoints)
	{
    	istensor ep;

    	// computes strain minus the initial strain, ie, the one that
    	// creates stresses
    	ev->computeStrain(dofset::tn1, ep);
    	strain += ep;
    	ev->updateCurrentState(when);
        
        
        istensor sp;
        ev->theMaterialPoint->stress(sp);
        stress += sp;
	}

    strain = strain * (1.0/ theEvaluationPoints.size());
    theState->exx = strain(0,0);
    theState->eyy = strain(1,1);
	theState->ezz = strain(2,2);
	theState->exy = strain(0,1);
	theState->eyz = strain(1,2);
	theState->exz = strain(0,2);


    stress = stress * (1.0/ theEvaluationPoints.size());
    theState->sxx = stress(0,0);
    theState->syy = stress(1,1);
	theState->szz = stress(2,2);
	theState->sxy = stress(0,1);
	theState->syz = stress(1,2);
	theState->sxz = stress(0,2);
}






