/*
 * Cluster.h
 *
 *  Created on: May 30, 2011
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

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include "Defect.h"
#include "ClusterParam.h"
#include "kernel/ParticleType.h"
#include "kernel/Coordinates.h"
#include "EDType.h"
#include <vector>
#include <map>
#include <string>
#include "Particle.h"

namespace Kernel { class Domain; class SubDomain; }

namespace OKMC {

class MobileParticle;

class Cluster: public OKMC::Defect {
public:
	Cluster(Kernel::SubDomain *, Kernel::Domain *, unsigned type, MobileParticle *, MobileParticle *);
	Cluster(Kernel::SubDomain *, Kernel::Domain *, unsigned type, const Kernel::Coordinates &);
	Cluster(std::istream &);

	Cluster * growDefect(Kernel::SubDomain *, const Kernel::ID &theMap);

	virtual ~Cluster();

	virtual void perform (Kernel::SubDomain *pSub, unsigned eventType);
	virtual void setIndex(unsigned eventType, int idx) { _idx[eventType] = idx; }
	virtual unsigned getNEvents() const { return Kernel::MAX_PARTICLES+4; } //(emission, migration, 2 transformations and recombination)
	virtual int  getIndex(unsigned eventType) const { return _idx[eventType]; }
	virtual float getRate(unsigned eventType, float kT) const;
	virtual E_TYPE getEType() const { return CLUSTER; }
	unsigned       getEDType() const { return _edtype; }

	virtual Defect * interact(Kernel::SubDomain *, Defect *dyn, Particle *sta, std::vector<Particle *> &);
	virtual bool canInteract (Kernel::SubDomain *, const Particle *dyn, const Particle *sta) const;
	virtual Kernel::ID getID() const { return _map; }
	virtual std::vector<Particle *> getParticles() { return _particles; }
	virtual Kernel::MeshElement * getElement() const { return _pElement; }
	virtual void deletePart(Kernel::SubDomain *, std::vector<Particle *> &);

	const EDType::CLType * getMyParam();
	unsigned getHash() const { return _hash; }
	virtual unsigned getState() const {return 0;}
	double getRadius() const;

	static double emissionRate     (const Kernel::MeshElement *, unsigned edtype, const EDType::CLType *, const Kernel::ID &, unsigned pt, double kT);
	static double recombinationRate(const Kernel::MeshElement *, unsigned edtype, const EDType::CLType *, const Kernel::ID &, double kT);

	static void addMap(Kernel::ID &first, const Kernel::ID &second, bool ivmodel);
	static void addTo(Kernel::P_TYPE pt, Kernel::ID &, bool ivmodel);
	static bool recombination(Kernel::ID &);
	static bool emitFrom(Kernel::P_TYPE pt, Kernel::ID &); // false if it cannot emit the particle V4 -> He + ???
	static void collapseMap(Kernel::ID &); //remove null entries.

	virtual void restart(std::ostream &) const;

private:
	int _idx[Kernel::MAX_PARTICLES+4];

	std::vector<Particle *> _particles;
	Kernel::ID _map;
	const EDType::CLType * _myParams;
	unsigned _hash;
	unsigned _edtype;
	Kernel::Coordinates _axes[3]; //0&1 parallel, 2 is perpendicular to the plane
	Kernel::Coordinates _center;
	float _radius;
	Kernel::MeshElement *_pElement;

	Kernel::Coordinates addParticle(Kernel::SubDomain *pSub, Particle *p); //returns the element where the
	//particle was tried to be added
	bool fuzzModel(Kernel::SubDomain *, Kernel::P_TYPE, const Kernel::Coordinates &, bool bKillIfInterface); //applies the fuzz model. Element is the one
	 // returned by addParticle. True if cluster was removed. It adds one particle
	void removeParticle(Kernel::P_TYPE fam);
	void emit(Kernel::SubDomain *pSub, unsigned);

	void setAxes(Kernel::SubDomain *);
	void migrate(Kernel::SubDomain *);
	void transform(Kernel::SubDomain *, unsigned);
	void recombine(Kernel::SubDomain *);
	Defect * checkSize(Kernel::SubDomain *pSub, bool bInsert, std::vector<Particle *> &); //insert in scheduler?
	void updateRadius(const Particle *part, int adding);
};

}

#endif /* CLUSTER_H_ */
