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
 *  volumeloading.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 9/9/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "volumeloading.h"

#include "loading.h"
#include "Analysis/assembler.h"
#include "Analysis/Propfactors/propfactor.h"
#include "Elements/evalspot.h"
#include "Model/model.h"
#include "Model/Node/node.h"
#include "Model/Parts/modelpart.h"
#include "Io/usercommand.h"
#include "Io/logger.h"

#include <boost/foreach.hpp>

#ifdef WITHTBB
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/parallel_for_each.h"
#endif

using namespace blue;



volumeLoading :: volumeLoading(const commandLine &cl) :
	loading(cl),
	theVolumeName(),
	theVolume(0)
{
	theType = LOAD_BODY;
	
	for (int j=1; j< cl.size(); j++)
	{
		const usercommand &uc = cl[j];
		if ( uc.keyword() == "volume" )	theVolumeName = uc.option();
	}
	
}




#ifdef WITHTBB

// data structure for tbb threaded version
class tbbAccumulateEnergySpots
{
public:
	energies my_energy;
	
	tbbAccumulateEnergySpots(vector<evalspot*>& ev, const double& when, const factorcombo* factor, const ivector& fU) :
        my_spots(ev), my_time(when), my_scaling(factor), my_fU(fU) {my_energy.setZero();}
    
	tbbAccumulateEnergySpots( tbbAccumulateEnergySpots& nn, tbb::split) :
        my_spots(nn.my_spots), my_time(nn.my_time), my_scaling(nn.my_scaling), my_fU(nn.my_fU) {my_energy.setZero();}
    
	void join(const tbbAccumulateEnergySpots& y)
    {
        my_energy.potential += y.my_energy.potential;
        my_energy.energy    += y.my_energy.energy;
    }
	
    
	void operator()( const tbb::blocked_range<size_t>& r)
	{
        const size_t re = r.end();
        
        for (int e = r.begin(); e != re; e++)
        {
            evalspot* ev = my_spots[e];
            ivector   coor = ev->coordinates();
            
			double pot = my_scaling->eval(my_time, coor) * my_fU.dot(coor) * ev->mass();
            my_energy.potential -= pot;
            my_energy.energy    -= pot;
		}
	}
	
private:
    std::vector<evalspot*>&         my_spots;
    const double                    my_time;
    const factorcombo*              my_scaling;
    const ivector&                  my_fU;
};




class tbbAccumulateEnergyElements
{
public:
	energies tenergy;
	
	tbbAccumulateEnergyElements( vector<element*>& ev, const dofset::evaluation_time& when) :
    my_elements(ev), evtime(when)
    {
        tenergy.setZero();
    }
    
	tbbAccumulateEnergyElements( tbbAccumulateEnergyElements& nn, tbb::split) : my_elements(nn.my_elements), evtime(nn.evtime)
    {
        tenergy.setZero();
    }
    
	void join(const tbbAccumulateEnergyElements& y)
    {
        tenergy += y.tenergy;
    }
	
	void operator()( const tbb::blocked_range<size_t>& r)
	{
        const size_t re = r.end();
        
        for (size_t e= r.begin(); e != re; e++)
        {
			energies elenergy;
			my_elements[e]->integrateEnergy(elenergy, evtime);
			tenergy += elenergy;
		}
	}
	
private:
    std::vector<element*>&          my_elements;
    const dofset::evaluation_time   evtime;
};

#endif




void volumeLoading :: accumulateEnergy(energies& ener, const double time)
{
	if (factor == 0) factor = &(factorcombo::getPropfactorCombinationFromLabel(getPropfactorLabel()));
    
    
#ifdef WITHTBB
    
	tbbAccumulateEnergySpots tbbenergy(theVolume->evalspots, time, factor, fU);
	tbb::parallel_reduce(tbb::blocked_range<size_t>(0, theVolume->evalspots.size()), tbbenergy);
    ener += tbbenergy.my_energy;
	
/*III
    tbbComputeEnergyElements tbbenergyE(bb->elements, when);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, elements.size()), tbbenergyE);
    energy += tbbenergyE.tenergy;
*/
    
#else
    BOOST_FOREACH(evalspot* ev, theVolume->evalspots)
    {
        double pot = factor->eval(time, ev->coordinates()) * fU.dot( ev->coordinates() ) * ev->mass();
        ener.potential -= pot;
        ener.energy    -= pot;
    }
#endif
}




void volumeLoading :: assembleLoading(assembler& asb, const double time, const double f)
{
	if (factor == 0) factor = &(factorcombo::getPropfactorCombinationFromLabel(getPropfactorLabel()));
 
    std::vector<element*>::iterator eliter = theVolume->elements.begin();
    while (eliter != theVolume->elements.end())
    {
        
        std::set<evalspot*> ev = (*eliter)->getEvalspots();
        
        std::set<evalspot*>::iterator eviter = ev.begin();
        while (eviter != ev.end())
        {
            std::vector<shapefun>& theShp  = (*eviter)->getShapefunctions();
            ivector force = fU * factor->eval(time, (*eviter)->coordinates() ) * (*eviter)->mass();
            
            for (size_t a=0; a< theShp.size(); a++)
            {
                Unode* nd = dynamic_cast<Unode*>(theShp[a].parentNode);
                asb.assemble( force * theShp[a].value , nd->getUDS() );
            }
            
            ++eviter;
        }
        
        ++eliter;
    }
}




void volumeLoading :: initialize(model& theModel)
{
	loading::initialize(theModel);
	theVolume = &(theModel.getPart(theVolumeName));
}




void volumeLoading :: print(std::ostream& of) const
{
    of  << "\n\n Volume loading on volume " << theVolumeName
        << "\n     f       : " << fU
        << "\n     Scaling : " << getPropfactorLabel();
}




void   volumeLoading :: transferLoadsToNodes(const double time)
{
	// pointloads are not initialized at the beginning and the pointer to the factor
	// is left to null. The first time we use them, we initilize the pointer, for faster
	// access
	if (factor == 0) factor = &(factorcombo::getPropfactorCombinationFromLabel(getPropfactorLabel()));
    
    std::vector<evalspot*>::iterator eviter = theVolume->evalspots.begin();
   while (eviter != theVolume->evalspots.end())
    {
        std::vector<shapefun>& theShp  = (*eviter)->getShapefunctions();
        
        std::vector<shapefun>::iterator shiter = theShp.begin();
        while ( shiter != theShp.end() )
        {
            (*shiter).parentNode->assembleForce( factor->eval(time, (*eviter)->coordinates() ) * (*eviter)->mass() * fU * (*shiter).value   );
            ++shiter;
        }
        ++eviter;
    }

    
    std::vector<element*>::iterator eliter = theVolume->elements.begin();
    while (eliter != theVolume->elements.end())
    {
       std::set<evalspot*> ev = (*eliter)->getEvalspots();
        std::set<evalspot*>::iterator eviter = ev.begin();
        while (eviter != ev.end())
        {
            std::vector<shapefun>& theShp  = (*eviter)->getShapefunctions();
            
            std::vector<shapefun>::iterator shiter = theShp.begin();
            while ( shiter != theShp.end() )
            {
                (*shiter).parentNode->assembleForce( factor->eval(time, (*eviter)->coordinates() ) * (*eviter)->mass() * fU * (*shiter).value   );

                ++shiter;
            }
            ++eviter;
        }
            
        ++eliter;
    }
}


