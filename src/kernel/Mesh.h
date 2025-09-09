/*
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

#ifndef KMCMESH_H
#define KMCMESH_H

#include <vector>
#include <unordered_map>
#include "Coordinates.h"

#include "MeshElement.h"

#include "kernel/ParticleType.h"

#include "electrostatics/MeshNode.h"
#include "electrostatics/ParticleToNodeHandler.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace LKMC { class LatticeSite; }
namespace OKMC { class Particle; class MobileParticle; }
namespace Domains { class MCClient; }

namespace Kernel {
class Domain;

class Mesh
{
public:
	struct LSNeiInfo
	{
		LSNeiInfo(LKMC::LatticeSite *p, float d2) : _pLS(p), _dist2(d2) {}
		LKMC::LatticeSite * _pLS;
		float _dist2;
		bool operator<(const LSNeiInfo &li) const { return _dist2 < li._dist2; }
	};
	struct ClusterInfo
	{
		OKMC::Particle * _pPart;
		Coordinates _coord;
		Coordinates _orig;
		MeshElement *_pEle;
	};
	typedef std::vector<LSNeiInfo> LSNeiList;
	typedef std::vector<MeshElement>::iterator iterator;
	
    Mesh(Domain *, const Coordinates &m, const Coordinates &M,
         const std::vector<float> &x, const std::vector<float> &y, const std::vector<float> &z,
	 const Domains::MCClient *);
    ~Mesh();

    void insert(LKMC::LatticeSite *, MeshElement *); //0 if element not known
    void remove(LKMC::LatticeSite *);
    void insert(OKMC::Particle *, MeshElement *);  // if new mesh element not known, put 0.
    void remove(OKMC::Particle *);

    void     getDomain(Coordinates &m, Coordinates &M) const;
    unsigned getIndexFromCoordinates(const Coordinates &) const;
    unsigned getIndexFromIndices(unsigned ix, unsigned iy, unsigned iz) const { return ix*_yzCells + iy*_zCells + iz; }
    void     getIndicesFromIndex(unsigned idx, unsigned &, unsigned &, unsigned &) const;
    MeshElement *getElement(unsigned i) { return &_elements[i]; }
    unsigned size()  { return _elements.size();  }
    iterator begin() { return _elements.begin(); }
    iterator end()   { return _elements.end();   }
    
    unsigned getxIndex(float) const;
    unsigned getyIndex(float) const;
    unsigned getzIndex(float) const;
    
    void getCenter(Coordinates &c, unsigned idx) const;
    void getCorners(unsigned, Coordinates &m, Coordinates &M) const;
    void getRandomPosition(unsigned, Coordinates &m, float, float, float) const;
    double getVolume(unsigned idx) const;
    const std::vector<float> & getLines(unsigned h) const
    		{ if(h==0) return _xlines; if(h==1) return _ylines; return _zlines; }
    void fillInterfaces(const Coordinates &, float dist, std::vector<OKMC::Interface *> &);

    bool fillLatticeNeighbors       (const Coordinates &, float dist, LSNeiList &); //LatticeAtoms
    void fillNeighborsAllMat        (const Coordinates &, float dist, std::vector<MeshElement *> &); //MeshElements
    void fillNeighborsTopology      (MeshElement *pEle, std::vector<const MeshElement *> &, std::vector<const MeshElement *> &, std::vector<const MeshElement *> &) const; //Edges, etc...
    void fillNeighborsOneMat        (const MeshElement *pEle, const Coordinates &, float dist, std::vector<MeshElement *> &elems);
    void fillNeighborsOneMat        (const MeshElement *pEle, const Coordinates &, float dist, std::vector<OKMC::Particle *> &parts);
    void fillAdjacentNeighbors      (const MeshElement *pEle, std::vector<MeshElement *> &all, M_TYPE mt, M_TYPE basicMt);
    void fillAdjacentNeighborsOneMat(const MeshElement *pEle, std::vector<MeshElement *> &); //adjacent Meshelements

    void setPeriodicRelative(const Coordinates &c, Coordinates &n) const; // n = n-c
    void setPeriodicRelative(const ublas::vector<float> &, ublas::vector<float> &) const;
    void setPeriodicRelative(const ublas::vector<double> &, ublas::vector<double> &) const;
    bool isInDomain(const Coordinates & c) const;
    void middle(Coordinates &c, const Coordinates &c1, const Coordinates &c2) const;
    LKMC::LatticeSite * findLS(const Coordinates &c, float dist=0.01);
    enum JUMP_ACTIONS { JUMP_OK, JUMP_REJECTED, JUMP_ANNIHILATED, JUMP_INTERFACE };
    JUMP_ACTIONS jumpPosition(Coordinates &c, const Coordinates &lambda, MeshElement *&pEle,
    		const std::pair<P_TYPE , unsigned> *jumpType, Coordinates &orig, unsigned lhf);
    JUMP_ACTIONS checkMove   (MeshElement * &pEle, Coordinates &to);
    void getInteractions       (SubDomain *, const OKMC::Particle *, std::vector<OKMC::Particle *> &parts);

    void resetDiffusivity();
    void changeMaterial(SubDomain *, unsigned idx, M_TYPE to);
    M_TYPE getMaterial(Coordinates const& aWhere) const;

    void print() const;
    void printDomains() const;

    Electrostatics::MeshNode              ***getNodes() const { return _pNodes; }
    Electrostatics::ParticleToNodeHandler   *getParticleToNodeHandler() const { return _pPNH; }

    unsigned getnx() const { return _xlines.size(); }
    unsigned getny() const { return _ylines.size(); }
    unsigned getnz() const { return _zlines.size(); }

    std::vector<float> const& getLinesX() const { return _xlines; }
    std::vector<float> const& getLinesY() const { return _ylines; }
    std::vector<float> const& getLinesZ() const { return _zlines; }

    bool getPeriodicX() const { return _periodicX; }
    bool getPeriodicY() const { return _periodicY; }
    bool getPeriodicZ() const { return _periodicZ; }

    void getMidX(unsigned, unsigned, unsigned, double &, double &);
    void getMidY(unsigned, unsigned, unsigned, double &, double &);
    void getMidZ(unsigned, unsigned, unsigned, double &, double &);

    void updatePotential(std::set<Electrostatics::MeshNode *> &);
    void updatePotential();

    bool isDirichletX() const;
    bool isDirichletY() const;
    bool isDirichletZ() const;

    void getElementsFromNode(Electrostatics::MeshNode *, std::set<MeshElement *> &) const;
    void getNodesFromElement(MeshElement *, std::set<Electrostatics::MeshNode*> &) const;
    Electrostatics::MeshNode * getFirstNodeFromElement(MeshElement *) const;
    
    unsigned longHopFactor(unsigned idx) const;
    void     trackLHF     (unsigned idx, int plus_minus);

    void changeNumberOfCells(Kernel::SubDomain *, MeshElement *, int howMany); //one cell is "expanding" or contracting. Try to create
    // and extension/hole at the surface.

    // Selfdiffusion for a AB alloy. A jump or B jump
    void jumpBorA(MeshElement * from, MeshElement * to, double prob);
    
    void gatherVtkMaterialData(float const aVtkCellDecrement, VtkMaterialData &aData) const;

private:
    const Coordinates _min, _max;
    Domain * _pDomain;

    std::vector<MeshElement> _elements;

    Electrostatics::MeshNode      ***_pNodes;
    Electrostatics::ParticleToNodeHandler *_pPNH;

    std::map<Electrostatics::MeshNode *, std::set<MeshElement *> > _MN2ME;
    std::map<MeshElement *, std::set<Electrostatics::MeshNode *> > _ME2MN;

    std::vector<float> _xlines;
    std::vector<float> _ylines;
    std::vector<float> _zlines;
    
    std::vector<unsigned> _xParticles;
    std::vector<unsigned> _yParticles;
    std::unordered_map<unsigned, double> _longHopFactor[MAX_MATERIALS];

    unsigned _xCells; //sizeZ
    unsigned _yCells; //sizeZ
    unsigned _zCells; //sizeZ
    unsigned _yzCells; //size Y * Z
    bool _periodicX, _periodicY, _periodicZ;
    double _xsize, _ysize, _zsize;
    double _xModOffset, _yModOffset, _zModOffset;

	std::string _poissonX, _poissonY, _poissonZ;
	double      _potentialX, _potentialY, _potentialZ;

    void buildNodes();
    void initNodes();
    void buildNeighbors();

    const MeshElement * step(unsigned ix, unsigned iy, unsigned iz, int xStep, int yStep, int zStep) const;

    void check() const;
    //transfer the list of particle in the set between materials.
    void transfer(Kernel::SubDomain *, std::vector<OKMC::Particle *> &pPart, MeshElement *pME, M_TYPE from, M_TYPE to);
};

}
#endif
