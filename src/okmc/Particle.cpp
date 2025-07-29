/*
 * Particle.cpp
 *
 *  Created on: Feb 24, 2011
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

#include "Particle.h"
#include "kernel/Constants.h"
#include "kernel/SubDomain.h"
#include "kernel/Domain.h"
#include "kernel/Mesh.h"
#include "kernel/MeshElement.h"
#include "kernel/MeshParam.h"
#include "kernel/RateManager.h"
#include "io/ParameterManager.h"
#include "okmc/AlloyParam.h"
#include "okmc/MobileParticle.h"
#include "okmc/MobileParticleParam.h"


namespace OKMC {

Particle::Particle(Kernel::P_TYPE type, const Kernel::Coordinates &c, Defect *pD, const Kernel::Coordinates &o) : _coord(c), _orig(o), _pDefect(pD), _ptype(type)
{
	_next  = _prev = 0;
	_pElement = 0;
}

Particle::Particle(std::istream &is, Defect *pD) : _pDefect(pD)
{
	int pt;
	is >> _coord >> _orig >> pt;
	_ptype = pt;
	_next  = _prev = 0;
	_pElement = 0;
}

void Particle::restart(std::ostream &os) const
{
	os << _coord << " " << _orig << " " << int(_ptype) << " ";
}

void Particle::selfdiffusion(Kernel::SubDomain *pSub, Kernel::MeshElement *pEle)
{
    Kernel::M_TYPE mt = _pElement->getMaterial();
    char iorv = Domains::global()->PM()->getIorV(mt, _ptype);
    
    if (iorv < 2 && Domains::global()->PM()->getMaterial(mt)._selfDiffusion)
    {
        Kernel::Coordinates from, to;
        Kernel::Domain * pDomain = _pElement->getDomain();
        Kernel::Mesh * pMesh = pDomain->_pMesh;
        pMesh->getCenter(from, _pElement->getIndex());
        pMesh->getCenter(to, pEle->getIndex());
        pMesh->setPeriodicRelative(from, to);
        
        IO::Arrhenius alpha_x1 = pDomain->_pAlPar->_alpha[mt][int(iorv)](_pElement);
        IO::Arrhenius alpha_x2 = pDomain->_pAlPar->_alpha[mt][int(iorv)](pEle);
        float T = pDomain->_pRM->getT();
        
        //			 *** SMOOTHING FOR ENERGIES ***
        std::vector<const Kernel::MeshElement *> sides1, edges1, corners1;
        std::vector<const Kernel::MeshElement *> sides2, edges2, corners2;
        pMesh->fillNeighborsTopology(_pElement, sides1, edges1, corners1);
        pMesh->fillNeighborsTopology(pEle, sides2, edges2, corners2);
        
        double wSides1 = pDomain->_pAlPar->_smoothing[_pElement->getMaterial()][0];
        double wEdges1 = pDomain->_pAlPar->_smoothing[_pElement->getMaterial()][1];
        double wCorners1 = pDomain->_pAlPar->_smoothing[_pElement->getMaterial()][2];
        double wCentral1 = 1 - 6 * wSides1 - 12 * wEdges1 - 8 * wCorners1;
        // Mirror boundary conditions
        wCentral1 += (6 - sides1.size()) * wSides1  + (12 - edges1.size()) * wEdges1 + (8 - corners1.size()) * wCorners1;
        
        double wSides2 = pDomain->_pAlPar->_smoothing[pEle->getMaterial()][0];
        double wEdges2 = pDomain->_pAlPar->_smoothing[pEle->getMaterial()][1];
        double wCorners2 = pDomain->_pAlPar->_smoothing[pEle->getMaterial()][2];
        double wCentral2 = 1 - 6 * wSides2 - 12 * wEdges2 - 8 * wCorners2;
        wCentral2 += (6 - sides2.size()) * wSides2  + (12 - edges2.size()) * wEdges2 + (8 - corners2.size()) * wCorners2;
        
        double dE1 = wCentral1 * pDomain->_pAlPar->getDerMixingEnergy(mt, _pElement->getEffectiveAlloyFraction(), T);
        double dE2 = wCentral2 * pDomain->_pAlPar->getDerMixingEnergy(mt, pEle->getEffectiveAlloyFraction(), T);
        // Computing smoothing in energies
        for(std::vector<const Kernel::MeshElement *>::iterator it = sides1.begin(); it != sides1.end(); ++it)
            dE1 += wSides1 * pDomain->_pAlPar->getDerMixingEnergy(mt, (*it)->getEffectiveAlloyFraction(), T);
        for(std::vector<const Kernel::MeshElement *>::iterator it = edges1.begin(); it != edges1.end(); ++it)
            dE1 += wEdges1 * pDomain->_pAlPar->getDerMixingEnergy(mt, (*it)->getEffectiveAlloyFraction(), T);
        for(std::vector<const Kernel::MeshElement *>::iterator it = corners1.begin(); it != corners1.end(); ++it)
            dE1 += wCorners1 * pDomain->_pAlPar->getDerMixingEnergy(mt, (*it)->getEffectiveAlloyFraction(), T);
        
        for(std::vector<const Kernel::MeshElement *>::iterator it = sides2.begin(); it != sides2.end(); ++it)
            dE2 += wSides2 * pDomain->_pAlPar->getDerMixingEnergy(mt, (*it)->getEffectiveAlloyFraction(), T);
        for(std::vector<const Kernel::MeshElement *>::iterator it = edges2.begin(); it != edges2.end(); ++it)
            dE2 += wEdges2 * pDomain->_pAlPar->getDerMixingEnergy(mt, (*it)->getEffectiveAlloyFraction(), T);
        for(std::vector<const Kernel::MeshElement *>::iterator it = corners2.begin(); it != corners2.end(); ++it)
            dE2 += wCorners2 * pDomain->_pAlPar->getDerMixingEnergy(mt, (*it)->getEffectiveAlloyFraction(), T);
        
        float kT = pDomain->_pRM->getkT();
        // a' from cell 1 to cell 2
        double alpha_prime12 = std::sqrt(alpha_x1.getRate(kT) * alpha_x2.getRate(kT)) *
        std::exp((dE1 - dE2) * .5 / kT);
        // a' from cell 2 to cell 1
        double alpha_prime21 = std::sqrt(alpha_x1.getRate(kT) * alpha_x2.getRate(kT)) *
        std::exp((dE2 - dE1) * .5 / kT);
        
        double x_prime12, x_prime21;  // x' from cell 1 to cell 2 and x' from cell 2 to cell 1
        
        x_prime12 = alpha_prime12 * _pElement->getBAtoms() / (_pElement->getAAtoms() + alpha_prime12 * _pElement->getBAtoms());
        x_prime21 = alpha_prime21 * pEle->getBAtoms() / (pEle->getAAtoms() + alpha_prime21 * pEle->getBAtoms());

        Kernel::Coordinates cFrom, cTo;
        pDomain->_pMesh->getCenter(cFrom, _pElement->getIndex());
        pDomain->_pMesh->getCenter(cTo, pEle->getIndex());

        // Geometrical factor for computing the atom jump probabilities
        double g = pDomain->_pAlPar->_corrFactor[mt][int(iorv)] *
                    	                pMesh->longHopFactor(_pElement->getIndex()) *
                    	                pDomain->_pMePar->_lambda[mt] / to.abs();
        if(iorv == Kernel::IV_I)             // Interstitial from 1 to 2
        {
        	double probB_12;
            if(!_pElement->getBAtoms())      // A from 1 to 2
            	probB_12 = 0.;
            else if(!_pElement->getAAtoms()) // B from 1 to 2
            	probB_12 = 1.;
            else if(!pEle->getBAtoms())      // Probability of a B atom jumping from 1 to 2 when 2 has no B atoms
            	probB_12 = x_prime12 * g;
            else if(!pEle->getAAtoms())    	 // Probability of a B atom jumping from 1 to 2 when 2 has no A atoms
            	probB_12 = 1 - (x_prime12 - 1) * g;
            else                             // Probability of a B atom jumping from 1 to 2
                probB_12 = x_prime12 * .5 * (1 + g) +
                x_prime21 * .5 * (1 - g);
            pDomain->_pMesh->jumpBorA(_pElement, pEle, probB_12);
        }
        else if (iorv == Kernel::IV_V)       // Vacancy from 1 to 2
        {
        	double probB_21;
            if(!pEle->getBAtoms())           // A from 2 to 1
            	probB_21 = 0.;
            else if(!pEle->getAAtoms())      // B from 2 to 1
            	probB_21 = 1.;
            else if(!_pElement->getBAtoms()) // Probability of a B atom jumping from 2 to 1 when 1 has no B atoms
            	probB_21 = x_prime21 * g;
            else if(!_pElement->getAAtoms()) // Probability of a B atom jumping from 2 to 1 when 1 has no A atoms
            	probB_21 = 1 - (1 - x_prime21) * g;
            else                             // Probability of a B atom jumping from 2 to 1
                probB_21 = x_prime21 * .5 * (1 + g) +
                x_prime12 * .5 * (1 - g);
            pDomain->_pMesh->jumpBorA(pEle, _pElement, probB_21);
        }
        else
        	ERRORMSG("[Selfdiffusion] MobileParticle different from I or V trying to move lattice atoms");
    }
}

}



