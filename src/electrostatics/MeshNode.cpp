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

#include "kernel/Domain.h"

#include "io/Diagnostic.h"

#include "MeshNode.h"

using namespace boost::numeric;
using namespace Kernel;
using namespace Electrostatics;

MeshNode::MeshNode() {
	_index     = 0;

	_ix        = 0;
	_iy        = 0;
	_iz        = 0;

	_xm        = 0;
	_xp        = 0;
	_ym        = 0;
	_yp        = 0;
	_zm        = 0;
	_zp        = 0;

	_volume    = 0;

	_V         = 0.;
	_Ec        = 0.;
	_Ev        = 0.;
	_Eg        = 0.;

	_eDensity  = 0;
	_hDensity  = 0.;

	_N_okmc    = 0.;  // FIXME

	_E         = ublas::zero_vector<double>(3);

	_stress    = ublas::zero_vector<double>(6);
	_strain    = ublas::zero_vector<double>(6);

	// _potential = 0;
	_active    = false;

	_Dm = _Dp = _D0 = 0;
}

MeshNode::~MeshNode() {

}

double MeshNode::getWeight(LKMC::LatticeAtom *pLA) const {
	std::map<LKMC::LatticeAtom *, double>::const_iterator it;

	it = _mLA.find(pLA);
	if (it != _mLA.end())
		return it->second;
	else
		return 0.;
}

double MeshNode::getWeight(OKMC::Particle *pPart) const {
	std::map<OKMC::Particle *, double>::const_iterator it;

	it = _mPart.find(pPart);
	if (it != _mPart.end())
		return it->second;
	else
		return 0.;

}

void MeshNode::insert(OKMC::Particle *pPart, double w) {
	std::map<OKMC::Particle *, double>::iterator it;

	it = _mPart.find(pPart);
	if (it != _mPart.end())
		it->second += w;
	else
		_mPart.insert(std::pair<OKMC::Particle *, double>(pPart, w));
}

void MeshNode::insert(LKMC::LatticeAtom *pLA, double w) {
	std::map<LKMC::LatticeAtom *, double>::iterator it;

	it = _mLA.find(pLA);
	if (it != _mLA.end())
		it->second += w;
	else
		_mLA.insert(std::pair<LKMC::LatticeAtom *, double>(pLA, w));
}

void MeshNode::remove(OKMC::Particle *pPart) {
	_mPart.erase(pPart);
}

void MeshNode::remove(LKMC::LatticeAtom *pLA) {
	_mLA.erase(pLA);
}
