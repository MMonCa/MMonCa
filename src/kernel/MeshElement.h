/*
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

#ifndef KMCMESHELEMENT_H
#define KMCMESHELEMENT_H

#include "Material.h"
#include "kernel/ParticleType.h"
#include "lkmc/LKMCMode.h"
#include "Coordinates.h"
#include "io/Diagnostic.h"
#include <vector>
#include <map>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric;

namespace LKMC { class LatticeSite; }
namespace OKMC { class Particle; class Interface; }
namespace Kernel
{
class Mesh;
class Domain;
class SubDomain;
using VtkMaterialData = std::map<Coordinates, M_TYPE, VtkOrderLess<Coordinates::type>>;

class MeshElement
{
public:
    MeshElement(Domain *);
    ~MeshElement();
    
    unsigned getIndex() const { return _index; }
    M_TYPE getMaterial() const { return _mat; }
    M_TYPE getToMat(LKMC::LKMCMode);
    OKMC::Particle * getFirstPart() const { return _firstPart; }
    LKMC::LatticeSite * getFirstLS() const { return _firstLS; }
    std::vector<OKMC::Interface *> & getInterfaces() { return _interfaces; }
    const std::vector<OKMC::Interface *> & getInterfaces() const { return _interfaces; }
    unsigned getCrystallineLA() const { return _crystallineLA; }
    unsigned getNonCrystallineLA() const { return _nonCrystallineLA; }
    void incLA(bool b) { if(b) _crystallineLA++; else _nonCrystallineLA++; }
    void decLA(bool b) { if(b) _crystallineLA--; else _nonCrystallineLA--; }
    void updateLAStatus  (SubDomain *, LKMC::LKMCMode);   //false -> true
    void updateUnLAStatus(SubDomain *, LKMC::LKMCMode); //true -> false
    
    // Binary alloy AB    
    unsigned getAtoms() const { return _AAtoms + _BAtoms; };
    unsigned getAAtoms() const { return _AAtoms; }
    unsigned getBAtoms() const { return _BAtoms; }
    void incAAtoms() { _AAtoms++; }
    void decAAtoms();
    void incBAtoms() { _BAtoms++; }
    void decBAtoms();
    double getAlloyFraction() const; // Return the fraction of B atoms in the cell
    double getBasisFraction() const; // Return the fraction of A atoms in the cell
    double getEffectiveAlloyFraction() const; // Return the effective fraction (smoothing) of B atoms in the cell

    void insert(OKMC::Interface *);
    void remove(OKMC::Interface *);

    Domain * getDomain() const { return _pDomain; }
    char getSubDomainIdx() const { return _subDomainIdx; }
    char getLevel() const { return _subDomainLevel; }

    const std::map<std::pair<P_TYPE , unsigned >, unsigned > & getHops() const
		{ return _hops; }

    double & strain_xx() { return _strain(0); }
    double   strain_xx() const { return _strain(0); }
    double & strain_yy() { return _strain(1); }
    double   strain_yy() const { return _strain(1); }
    double & strain_zz() { return _strain(2); }
    double   strain_zz() const { return _strain(2); }
    
    double & strain_xy() { return _strain(3); }
    double   strain_xy() const { return _strain(3); }
    double & strain_xz() { return _strain(4); }
    double   strain_xz() const { return _strain(4); }
    double & strain_yz() { return _strain(5); }
    double   strain_yz() const { return _strain(5); }
    
    double & stress_xx() { return _stress(0); }
    double   stress_xx() const { return _stress(0); }
    double & stress_yy() { return _stress(1); }
    double   stress_yy() const { return _stress(1); }
    double & stress_zz() { return _stress(2); }
    double   stress_zz() const { return _stress(2); }
    
    double & stress_xy() { return _stress(3); }
    double   stress_xy() const { return _stress(3); }
    double & stress_xz() { return _stress(4); }
    double   stress_xz() const { return _stress(4); }
    double & stress_yz() { return _stress(5); }
    double   stress_yz() const { return _stress(5); }
    
    double   strain(unsigned i) const { return _strain(i); }
    ublas::vector<double> & strain()  { return _strain; }
    
    double   stress(unsigned i) const { return _stress(i); }
    ublas::vector<double> & stress()  { return _stress; }
    
    double & electrostaticPotential() { return _V; }
    double   electrostaticPotential() const { return _V; }
    double & bandGap() { return _Eg; }
    double   bandGap() const { return _Eg; }
    
    unsigned getHops(P_TYPE pt, unsigned state) const;
    void resetHops() { _hops.clear(); }
    
    unsigned getNumberOfParticles() const { return _howManyParts; }
    Kernel::Coordinates getCoords(float, float, float) const;
    void getCorners(Kernel::Coordinates &m, Kernel::Coordinates &M) const;
    double getVolume() const;

    void amorphizeME(Kernel::SubDomain *pSub);
    void setAmorphParts(float nParts) { _amorphParts = nParts; }
    void incrAmorphParts() { _amorphParts++; }
    float getAmorphParts() const {return _amorphParts;}

    void restart(std::ostream &) const;
    void restart(std::istream &);

    void gatherVtkMaterialData(float const aVtkCellDecrement, VtkMaterialData &aData) const;
		bool isCloseToGas(Kernel::Coordinates const& aWhere, float const aMargin) const;

    int _ABalance; // This is the balance between A atoms emitted/absorbed
    int _BBalance; // This is the balance between B atoms emitted/absorbed

private:
	LKMC::LatticeSite *_firstLS; //8
	OKMC::Particle    *_firstPart; //8
	unsigned _howManyParts;
	std::vector<OKMC::Interface *> _interfaces; //24
	std::vector<const MeshElement *> _sides; //32
	std::vector<const MeshElement *> _edges; //96
	std::vector<const MeshElement *> _corners; //64

	M_TYPE _mat; //4
	Domain *_pDomain; //8
	unsigned _index;
	char _subDomainIdx;
	char _subDomainLevel;
	
	ublas::vector<double> _strain; //24
	ublas::vector<double> _stress; //24

	double _V;   //8
	double _Eg;  //8        
	
	unsigned _nonCrystallineLA; //8
	unsigned _crystallineLA; //8
    unsigned _AAtoms; //8
	unsigned _BAtoms; //8
	float _amorphParts;

	//                   pt   state -1    = 8 hops
	std::map<std::pair<P_TYPE, unsigned >, unsigned > _hops;  //48
        
	friend class Mesh;
};

}

#endif
