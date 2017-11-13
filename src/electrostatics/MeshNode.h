/* Author: Benoit Sklenard benoit.sklenard@cea.fr 
 * 
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

#ifndef KMCMESHNODE_H
#define KMCMESHNODE_H

// #include "okmc/Particle.h"
#include "kernel/Coordinates.h"

#include "kernel/ParticleType.h"
#include "lkmc/LatticeAtom.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace OKMC { class Particle; class Interface; }
namespace LKMC { class LatticeAtom; }

// class Kernel::Mesh;
// class Kernel::Domain;

namespace Electrostatics
{

class MeshNode
{
	long int    _index;
	unsigned    _ix;
	unsigned    _iy;
	unsigned    _iz;

        Kernel::Coordinates _coord;

	boost::numeric::ublas::vector<double> _stress;
	boost::numeric::ublas::vector<double> _strain;

	boost::numeric::ublas::vector<double> _E; // Electric field
	double      _V;           // Electrostatic potential (eV)
	double      _Ec;          // CBM (eV)
	double      _Ev;          // VBM (eV)
	double      _Eg;          // Gap (usefull with a BGN model) (eV)
	double      _Eg0;         // Gap (at 0K) (eV)
	double      _eDensity;    // Electron density (/nm_3)
	double      _hDensity;    // Hole density (/nm⁻3)
	double      _Nc;          // Electron density of states (/nm3)
	double      _Nv;          // Holes density of states

	double      _N_okmc; // FIXME

	double      _xm; //nm
	double      _xp;
	double      _ym;
	double      _yp;
	double      _zm;
	double      _zp;
	double      _volume; //nm3
	bool        _active;

	double _Dm;
	double _Dp;
	double _D0;

	// add the information of active/non active nodes

	std::map<OKMC::Particle *,    double> _mPart;
	std::map<LKMC::LatticeAtom *, double> _mLA;

public:
	MeshNode();
    ~MeshNode();

    void insert(OKMC::Particle *, double);
    void remove(OKMC::Particle *);

    void insert(LKMC::LatticeAtom *, double);
    void remove(LKMC::LatticeAtom *);

    double getPotential() const          { return _V; }
    double getBandgap() const            { return _Eg; }
    double getBandgap0() const           { return _Eg0; }
    double getConductionBandedge() const { return _Ec; }
    double getValenceBandedge() const    { return _Ev; }
    double getElectronDensity() const    { return _eDensity; }
    double getHoleDensity() const        { return _hDensity; }
    double getVolume() const             { return _volume; }

    double getDopantDensity() const      { return _N_okmc; }

    void   getCoordinates(Kernel::Coordinates &coord) const { coord = _coord; }

    boost::numeric::ublas::vector<double> & getElectricField()  { return _E; }

    double getWeight(LKMC::LatticeAtom *) const;
    double getWeight(OKMC::Particle *) const;

    friend class Kernel::Mesh;
    friend class ParticleToNodeHandler;
    friend class LocalFermi;
    friend class Poisson;
};
}

#endif /* ! KMCMESHNODE_H */
