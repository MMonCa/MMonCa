/*
 * Author: Benoit Sklenard benoit.sklenard@cea.fr 
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

#ifndef PARTICLETONODEHANDLER_H
#define PARTICLETONODEHANDLER_H

#include <map>
#include <set>
#include <vector>

#include "kernel/Domain.h"
#include "kernel/Mesh.h"

#include "MeshNode.h"
#include "okmc/Particle.h"
#include "lkmc/LatticeAtom.h"

// namespace OKMC { class Particle; class MobileParticle; }
// namespace LKMC { class LatticeAtom; }

// class Kernel::Mesh;

namespace Electrostatics
{

class ParticleToNodeHandler {
	std::map<OKMC::Particle *,    std::set<MeshNode *> > _syncOKMC;
	std::map<LKMC::LatticeAtom *, std::set<MeshNode *> > _syncLKMC;

        Kernel::Mesh   *_pMesh;
        Kernel::Domain *_pDomain;

	inline double getOverlap(double, double);
	double addWeightNode(MeshNode *, OKMC::Particle *);    // DO NOT USE
	double addWeightNode(MeshNode *, LKMC::LatticeAtom *); // DO NOT USE

public:
	ParticleToNodeHandler(Kernel::Domain *, Kernel::Mesh *);
	~ParticleToNodeHandler();

	void remove(OKMC::Particle *);
	void insert(OKMC::Particle *);

	void remove(LKMC::LatticeAtom *);
	void insert(LKMC::LatticeAtom *);

	void getParticleNodes(OKMC::Particle *, std::set<MeshNode *> &);
	void getParticleNodes(LKMC::LatticeAtom *, std::set<MeshNode *> &);

	void print();
	void printLKMC();
};

}

#endif /* ! PARTICLETONODEHANDLER_H */
