/*
 * LatticeDiamond.cpp
 *
 *  Created on: Aug 17, 2012
 *      Author: ignacio.martin@imdea.org
 */

#include "LatticeDiamond.h"
#include "LatticeDiamondParam.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"

using std::vector;
using std::pair;
using std::set;
using Kernel::Coordinates;

namespace LKMC {

LatticeDiamond::LatticeDiamond(Kernel::Domain *p, const LatticeParam * plp, Kernel::M_TYPE mt)
: Lattice(p, plp, mt )
{
	_twinProb     = static_cast<const LatticeDiamondParam *>(plp)->_twinProb;
	_distortion   = static_cast<const LatticeDiamondParam *>(plp)->_growthDistortion;
	setSystems();
}

LatticeDiamond::~LatticeDiamond() {

}

/* fill ignores the lattice for alloys, but getNei? uses it */
void LatticeDiamond::fill(const Coordinates &mc, const Coordinates &Mc,
		vector<LatticeInformation> &pos, bool onlyFrame) const
{
	int m[3], M[3];
	getIndices(mc, Mc, m, M);
	double a=_latticeParameter[0];
	for(int i=m[0]; i<=M[0]; ++i)
		for(int j=m[1]; j<=M[1]; ++j)
			for(int k=m[2]; k<=M[2]; ++k)
			{
				tryInBox(pos, mc, Mc, 1, 1, toRotated(Coordinates(a*(i    ), a*(j    ), a*(k    ))));
				if(onlyFrame)
					continue;
				tryInBox(pos, mc, Mc, 1, 1, toRotated(Coordinates(a*(i+.50), a*(j+.50), a*(k    ))));
				tryInBox(pos, mc, Mc, 0, 0, toRotated(Coordinates(a*(i+.25), a*(j+.25), a*(k+.25))));
				tryInBox(pos, mc, Mc, 0, 0, toRotated(Coordinates(a*(i+.75), a*(j+.75), a*(k+.25))));
				tryInBox(pos, mc, Mc, 1, 1, toRotated(Coordinates(a*(i+.50), a*(j    ), a*(k+.50))));
				tryInBox(pos, mc, Mc, 1, 1, toRotated(Coordinates(a*(i    ), a*(j+.50), a*(k+.50))));
				tryInBox(pos, mc, Mc, 0, 0, toRotated(Coordinates(a*(i+.75), a*(j+.25), a*(k+.75))));
				tryInBox(pos, mc, Mc, 0, 0, toRotated(Coordinates(a*(i+.25), a*(j+.75), a*(k+.75))));
			}
}

/* Fill the system without lattice parameter */
void LatticeDiamond::setSystems()
{
	for(int s = 0; s < 2; ++s)
	{
		double l = _latticeParameter[s];
		_systems[s].push_back(vector<Coordinates>());
		_systems[s].back().push_back(toRotated(Coordinates(-.25, -.25, -.25)*l));
		_systems[s].back().push_back(toRotated(Coordinates(-.25,  .25,  .25)*l));
		_systems[s].back().push_back(toRotated(Coordinates( .25, -.25,  .25)*l));
		_systems[s].back().push_back(toRotated(Coordinates( .25,  .25, -.25)*l));

		for(int i=1; i<10; ++i)
		{
			_systems[s].push_back(vector<Coordinates>());
			for(unsigned j=0; j<4; ++j)
				if(i%2 == 1) //odd, center reflection of the previous one
					_systems[s].back().push_back(_systems[s][i-1][j]*-1.);
				else //even, create new twin defect as line reflection
				{
					int origin = i/2 -1; //this one does not change
					Coordinates a = _systems[s][0][origin]; //line
					Coordinates v = _systems[s][0][j]; //to reflect
					_systems[s].back().push_back(a*(a*v)*2./(a*a) - v);
				}
		}

		MEDMSG("Crystalline systems---");
		for(unsigned i=0; i<_systems[s].size(); ++i)
		{
			for(int j=0; j<4; ++j)
				MEDMSG(_systems[s][i][j]);
			MEDMSG("----");
		}
	}
}

// Given two relative coordinates (i.e., seen as if the atom is in 0,0,0) for two bonds, it returns the coordinates of the other
//neighbors/bonds (also seen from the atom).
//alloy = 0 or 1
const float distRef2 = 0.15*0.15;

unsigned LatticeDiamond::getNeis2(Kernel::SubDomain *pSub, int alloy, const Coordinates &one, const Coordinates &two, vector<pair<Coordinates, unsigned> > &all) const
{
	assert(sqrt(one*one) <= .3 && sqrt(two*two) <= .3);
	for(unsigned i=0; i<_systems[alloy].size(); ++i)
	{
		unsigned howMany = 0;
		//compute the distances. How many are smaller than distRef2?
		for(int j=0; j<4; ++j)
		{
			const Coordinates &ref = _systems[alloy][i][j];
			float d2 = (ref._x - one._x)*(ref._x - one._x) + (ref._y - one._y)*(ref._y - one._y) + (ref._z - one._z)*(ref._z - one._z);
			if(d2 < distRef2)
			{
				howMany++;
				continue;
			}
			d2 = (ref._x - two._x)*(ref._x - two._x) + (ref._y - two._y)*(ref._y - two._y) + (ref._z - two._z)*(ref._z - two._z);
			if(d2 < distRef2)
			{
				howMany++;
				continue;
			}
		}
		assert(howMany < 3);
		if(howMany == 2)
		{
			for(int j=0; j<4; ++j)
			{
				all.push_back(pair<Coordinates, unsigned>(_systems[alloy][i][j], i));
				if(_distortion)
					all.back().first += Coordinates(pSub->_rng.rand(), pSub->_rng.rand(), pSub->_rng.rand()) * _distortion;
			}
			return i;
		}
	}
	//not found (defect)
	return 15;
}

//given one single atomic bond returns all the possible atomic positions... including all the twins!!
void LatticeDiamond::getNeis1(Kernel::SubDomain *pSub, int alloy, const Kernel::Coordinates &one, unsigned ori,
				std::vector<std::pair<Kernel::Coordinates, unsigned> > &normal, std::vector<std::pair<Kernel::Coordinates, unsigned> > &twins) const
{
	assert(sqrt(one*one) <= .3);
	for(unsigned i=0; i<_systems[alloy].size(); ++i)
	{
		//compute the distances. How many are smaller than distRef2?
		for(unsigned j=0; j<4; ++j)
		{
			const Coordinates &ref = _systems[alloy][i][j];
			float d2 = (ref._x - one._x)*(ref._x - one._x) + (ref._y - one._y)*(ref._y - one._y) + (ref._z - one._z)*(ref._z - one._z);
			if(d2 < distRef2)
			{
				for(int k=0; k<4; ++k)
					if(i/2 == ori/2)
					{
						normal.push_back(pair<Coordinates, unsigned>(_systems[alloy][i][k], i));
						if(_distortion)
							normal.back().first += Coordinates(pSub->_rng.rand(), pSub->_rng.rand(), pSub->_rng.rand()) * _distortion;
					}
					else
					{
						twins.push_back(pair<Coordinates, unsigned>(_systems[alloy][i][k], i));
						if(_distortion)
							twins.back().first += Coordinates(pSub->_rng.rand(), pSub->_rng.rand(), pSub->_rng.rand()) * _distortion;
					}
				j=4;
			}
		}
	}
}

//one is the relative vector from amorphous atom 1, and two from 2. It tries to see which system can accomodate
//a correct atom between both of them.
//it does so by computing the systems of both atoms and seeing with ones are common.
//if several, then it picks one, that might be a twin with the "correct" probability.
unsigned LatticeDiamond::getNeis1(Kernel::SubDomain * pSubD, int alloy, const Coordinates &one, const Coordinates &two,
		vector<pair<Coordinates, unsigned> > &all, unsigned ori0, unsigned ori1) const
{
	int indexTwins=-1, indexOKs=-1; //systems
	for(unsigned i=0; i<_systems[alloy].size(); ++i)
	{
		bool isHere[2] = {false, false};
		for(int j=0; j<4; ++j)
		{
			const Coordinates &ref = _systems[alloy][i][j];
			float d2 = (ref._x - one._x)*(ref._x - one._x) + (ref._y - one._y)*(ref._y - one._y) + (ref._z - one._z)*(ref._z - one._z);
			if(d2 < distRef2)
				isHere[0] = true;
			d2 = (ref._x - two._x)*(ref._x - two._x) + (ref._y - two._y)*(ref._y - two._y) + (ref._z - two._z)*(ref._z - two._z);
			if(d2 < distRef2)
				isHere[1] = true;
		}
		if(isHere[0] && isHere[1])
		{
			if(i/2 == ori0/2 || i/2 == ori1/2)
				indexOKs = i;
			else
				indexTwins = i;
		}
	}
	unsigned idx = 15;
	if(indexTwins != -1 && indexOKs != -1)
		idx = (pSubD->_rng.rand() < _twinProb? indexTwins : indexOKs);
	else if(indexTwins != -1)
		idx = indexTwins;
	else if(indexOKs != -1)
		idx = indexOKs;
	if(idx != 15)
		for(int j=0; j<4; ++j)
		{
			all.push_back(pair<Coordinates, unsigned>(_systems[alloy][idx][j], idx));
			if(_distortion)
				all.back().first += Coordinates(pSubD->_rng.rand(), pSubD->_rng.rand(), pSubD->_rng.rand()) * _distortion;
		}
	return idx;
}

} /* namespace LKMC */
