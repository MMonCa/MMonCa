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

#include "Lattice.h"
#include "LatticeAtom.h"
#include "io/Diagnostic.h"
#include "kernel/Domain.h"
#include <algorithm>
#include <cmath>
#include <utility>

using namespace LKMC;
using Kernel::Coordinates;
using std::vector;
using std::pair;
using std::set;

Lattice::Lattice(Kernel::Domain *pD, const LatticeParam *p, Kernel::M_TYPE mt) :
	_pDomain(pD),_mt(mt)
{

	_latticeParameter[0] = p->_latticeParameter[0];
	_latticeParameter[1] = p->_latticeParameter[1];

	setMatrices(p->_mIndex[0], p->_mIndex[1], p->_mIndex[2],
		    p->_fIndex[0], p->_fIndex[1], p->_fIndex[2]);
}

void Lattice::getIndices(const Coordinates &mc, const Coordinates &Mc, int m[3], int M[3]) const
{
	int i[8], j[8], k[8], l=0;
	Coordinates c = Coordinates(mc._x, mc._y, mc._z);
	toLatticeUnits(i[l], j[l], k[l], c);
	l++;
	c = Coordinates(Mc._x, mc._y, mc._z);
	toLatticeUnits(i[l], j[l], k[l], c);
	l++;
	c = Coordinates(mc._x, Mc._y, mc._z);
	toLatticeUnits(i[l], j[l], k[l], c);
	l++;
	c = Coordinates(Mc._x, Mc._y, mc._z);
	toLatticeUnits(i[l], j[l], k[l], c);
	l++;
	c = Coordinates(mc._x, mc._y, Mc._z);
	toLatticeUnits(i[l], j[l], k[l], c);
	l++;
	c = Coordinates(Mc._x, mc._y, Mc._z);
	toLatticeUnits(i[l], j[l], k[l], c);
	l++;
	c = Coordinates(mc._x, Mc._y, Mc._z);
	toLatticeUnits(i[l], j[l], k[l], c);
	l++;
	c = Coordinates(Mc._x, Mc._y, Mc._z);
	toLatticeUnits(i[l], j[l], k[l], c);
	m[0] = i[0], m[1] = j[0], m[2] = k[0], 
	M[0] = i[0], M[1] = j[0], M[2] = k[0];
	for(l=0; l<8; ++l)
	{
		if(m[0] > i[l]) m[0] = i[l];
		if(M[0] < i[l]) M[0] = i[l];
		if(m[1] > j[l]) m[1] = j[l];
		if(M[1] < j[l]) M[1] = j[l];
		if(m[2] > k[l]) m[2] = k[l];
		if(M[2] < k[l]) M[2] = k[l];
	}
	for(l=0; l<3; ++l)
	{
		m[l]--;
		M[l]++;
	}
}

void Lattice::tryInBox(vector<LatticeInformation> &pos,
		const Coordinates &mc, const Coordinates &Mc, int orientation, unsigned type, const Coordinates &c) const
{
	if(c.isInto(mc, Mc))
	{
		LatticeInformation li;
		li._coords = c;
		li._orientation = orientation;
		li._type = type;
		pos.push_back(li);
	}
}


void Lattice::toLatticeUnits(int &mi, int &mj, int &mk, const Coordinates &c) const
{
	Coordinates original = fromRotated(c);
	mi = int(original._x/_latticeParameter[0]);
	mj = int(original._y/_latticeParameter[0]);
	mk = int(original._z/_latticeParameter[0]);
}

Coordinates Lattice::toRotated(const Coordinates &o) const
{
	Coordinates n;
	n._x = _matR[0][0]*o._x + _matR[0][1]*o._y + _matR[0][2]*o._z;
	n._y = _matR[1][0]*o._x + _matR[1][1]*o._y + _matR[1][2]*o._z;
	n._z = _matR[2][0]*o._x + _matR[2][1]*o._y + _matR[2][2]*o._z;
	return n;
}

Coordinates Lattice::fromRotated(const Coordinates &n) const
{
	Coordinates o;
	o._x = _invR[0][0]*n._x + _invR[0][1]*n._y + _invR[0][2]*n._z;
	o._y = _invR[1][0]*n._x + _invR[1][1]*n._y + _invR[1][2]*n._z;
	o._z = _invR[2][0]*n._x + _invR[2][1]*n._y + _invR[2][2]*n._z;
	return o;
}

void Lattice::getMatR(boost::numeric::ublas::matrix<double> &m) const {
	for (unsigned i = 0; i < 3; ++i)
		for (unsigned j = 0; j < 3; ++j)
			m(i, j) = _matR[i][j];
}

void Lattice::getInvR(boost::numeric::ublas::matrix<double> &m) const {
	for (unsigned i = 0; i < 3; ++i)
		for (unsigned j = 0; j < 3; ++j)
			m(i, j) = _invR[i][j];
}

void Lattice::setMatrices(float Mi, float Mj, float Mk, float Fi, float Fj, float Fk)
{

	// some checks.
    if(fabs(Mi*Fi + Mj*Fj + Mk*Fk) > 1e-5)
		ERRORMSG("Surface orientation not perpendicular to flat orientation");
    const double modM = std::sqrt(Mi*Mi + Mj*Mj + Mk*Mk);
    const double modF = std::sqrt(Fi*Fi + Fj*Fj + Fk*Fk);
    Mi /= modM; Mj /=modM; Mk /= modM;
    Fi /= modF; Fj /=modF; Fk /= modF;

	const double det = (Fj*Fj*Mk*Mk + Fi*Fi*Mk*Mk -2*Fj*Fk*Mj*Mk-2*Fi*Fk*Mi*Mk+
			Fk*Fk*Mj*Mj+Fi*Fi*Mj*Mj-2*Fi*Fj*Mi*Mj+Fk*Fk*Mi*Mi+Fj*Fj*Mi*Mi);

	_invR[0][0] = Mi;
	_invR[0][1] = Fi;
	_invR[0][2] = Mj*Fk - Fj*Mk;
	_invR[1][0] = Mj;
	_invR[1][1] = Fj;
	_invR[1][2] = Fi*Mk - Mi*Fk;
	_invR[2][0] = Mk;
	_invR[2][1] = Fk;
	_invR[2][2] = Mi*Fj - Fi*Mj;

	_matR[0][0] = (Fj*(Fj*Mi-Fi*Mj)-Fk*(Fi*Mk-Fk*Mi))/det;
	_matR[0][1] = (Fk*(Fk*Mj-Fj*Mk)-Fi*(Fj*Mi-Fi*Mj))/det;
	_matR[0][2] = (Fi*(Fi*Mk-Fk*Mi)-Fj*(Fk*Mj-Fj*Mk))/det;
	_matR[1][0] = (Mk*(Fi*Mk-Fk*Mi)-Mj*(Fj*Mi-Fi*Mj))/det;
	_matR[1][1] = (Mi*(Fj*Mi-Fi*Mj)-Mk*(Fk*Mj-Fj*Mk))/det;
	_matR[1][2] = (Mj*(Fk*Mj-Fj*Mk)-Mi*(Fi*Mk-Fk*Mi))/det;
	_matR[2][0] = (Fk*Mj-Fj*Mk)/det;
	_matR[2][1] = (Fi*Mk-Fk*Mi)/det;
	_matR[2][2] = (Fj*Mi-Fi*Mj)/det;

	HIGHMSG("--");
	HIGHMSG("[" << _invR[0][0] << " " << _invR[0][1] << " " << _invR[0][2] << "]");
	HIGHMSG("[" << _invR[1][0] << " " << _invR[1][1] << " " << _invR[1][2] << "]");
	HIGHMSG("[" << _invR[2][0] << " " << _invR[2][1] << " " << _invR[2][2] << "]");
	HIGHMSG("--");
	HIGHMSG("[" << _matR[0][0] << " " << _matR[0][1] << " " << _matR[0][2] << "]");
	HIGHMSG("[" << _matR[1][0] << " " << _matR[1][1] << " " << _matR[1][2] << "]");
	HIGHMSG("[" << _matR[2][0] << " " << _matR[2][1] << " " << _matR[2][2] << "]");
}
