/*
 * MobileParticle.h
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

#ifndef MOBILEPARTICLE_H_
#define MOBILEPARTICLE_H_

#include "Particle.h"
#include "Defect.h"
#include "domains/Global.h"
#include "io/ParameterManager.h"
#include "kernel/MeshElement.h"

namespace OKMC {

class MobileParticle: public Defect, public Particle
{
public:
	MobileParticle(Kernel::SubDomain *, Kernel::P_TYPE type, unsigned state, Kernel::Domain *p, const Kernel::Coordinates &c, const Kernel::Coordinates &o);
	MobileParticle(std::istream &);

	virtual ~MobileParticle();

	virtual void setIndex(unsigned eventType, int idx) { _idx[eventType] = idx; }
	virtual int  getIndex(unsigned eventType) const { return _idx[eventType]; }
	virtual unsigned getNEvents() const { return 6; } //migrate, br1, br2, FTI, FTV, states
	virtual float getRate(unsigned eventType, float kT) const;
	virtual Event::E_TYPE getEType() const { return MOBILEPARTICLE; }
	virtual void perform (Kernel::SubDomain *, unsigned eventType);
	virtual Kernel::MeshElement * getElement() const { return _pElement; }

	virtual Defect * interact(Kernel::SubDomain *, Defect *dyn, Particle *sta, std::vector<Particle *> &);
	virtual bool canInteract (Kernel::SubDomain *, const Particle *dyn, const Particle *sta) const;

	virtual Kernel::ID getID() const;
	virtual std::vector<Particle *> getParticles() { return std::vector<Particle *>(1,this); }
	void deletePart(Kernel::SubDomain *, std::vector<Particle *> &);
	virtual unsigned getState() const { return _state; }
        bool    fits(Kernel::MeshElement const * const pElem) const { return Domains::global()->PM()->getStates(pElem->getMaterial(), _ptype) > _state; }

	virtual void     restart(std::ostream &) const;

private:
	int _idx[6];
	void migrate(Kernel::SubDomain *);
	mutable unsigned _longHopFactor;
	Kernel::Coordinates _axes;
	void setAxes(Kernel::SubDomain *); //initialize migration axis.
	double breakUpClusterRate(float kT, int ev) const;
	void breakup(Kernel::SubDomain *, Kernel::P_POS);
	void emit(Kernel::SubDomain *, unsigned char iv);
	void updateState(Kernel::SubDomain *);
	Defect * internalInteract(Kernel::SubDomain *pSub, Defect *def, bool bForce, std::vector<Particle *> &);
};

}

#endif /* MOBILEPARTICLE_H_ */
