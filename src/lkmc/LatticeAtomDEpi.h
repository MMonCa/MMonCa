/*
 * LatticeAtomDEpi.h
 *
 *  Created on: Jul 15, 2014
 *      Author: ignacio.martin@imdea.org
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

#ifndef LATTICEATOMDEPI_H_
#define LATTICEATOMDEPI_H_

#include "LatticeAtomDiamond.h"
#include "EpiGasParam.h"

namespace LKMC {

class LatticeAtomDEpi : public LatticeAtomDiamond
{
public:
	LatticeAtomDEpi(LASTATE, Kernel::Domain *p, Kernel::M_TYPE mt, Kernel::P_TYPE type, const Kernel::Coordinates &c, unsigned orientation);
	LatticeAtomDEpi(std::istream &);

	virtual void perform (Kernel::SubDomain *, unsigned eventType);
	virtual float getRate(unsigned eventType, float kT) const;
	virtual void setIndex(unsigned eventType, int idx) { _idx[eventType] = idx; }
	virtual int  getIndex(unsigned eventType) const { return _idx[eventType]; }
	virtual unsigned getNEvents() const { return MAX_LKMC_EVENT; } //depositions + mig + etching + final adsorption + desorption
	//              4         4+12=16   16+12=28         Perfect, but due to periodic BC and twins it can be more.
	enum { FIRSTN = 8, SECONDN = 28, THIRDN = 48 };
	virtual unsigned first() const  { return FIRSTN; }
	virtual unsigned second() const { return SECONDN; }
	virtual unsigned third() const  { return THIRDN; }

	virtual std::string getClass() const { return "DEPI"; }
	virtual void restart(std::ostream &) const;

private:
	int _idx[MAX_LKMC_EVENT];
	const static unsigned MAX_ATOM_MIGS = 5;
	mutable unsigned _rateMigs;  //number of atoms in the cache
	mutable float    _rateMig[MAX_ATOM_MIGS];  //Mig freq. for each atom
	mutable LatticeAtomDEpi *_atomMig[MAX_ATOM_MIGS]; //position to go...
	Kernel::Coordinates _twinPos;
	bool _bTwin;

	double getMigRate(float formEner, float kT) const;
	void performEtching(Kernel::SubDomain *pSub, std::map<int, LatticeAtomDiamond *> &);
	void performPrecursor(Kernel::SubDomain *, Kernel::P_TYPE, std::map<int, LatticeAtomDiamond *> &);
	void performMig(Kernel::SubDomain *, std::map<int, LatticeAtomDiamond *> &);
	void performAdsorption(Kernel::SubDomain *, std::map<int, LatticeAtomDiamond *> &);
	void performDesorption(Kernel::SubDomain *, std::map<int, LatticeAtomDiamond *> &);
	//void performTwin(Kernel::SubDomain *, std::map<int, LatticeAtomDiamond *> &);
	void fillWithAtoms(Kernel::SubDomain *, LASTATE, std::vector<std::pair<Kernel::Coordinates, unsigned> > &neis,
			           std::vector<std::pair<Kernel::Coordinates, unsigned> > &twins, std::map<int, LatticeAtomDiamond *> &);
	void fillNeighborBoxes(Kernel::SubDomain *pSub);

	static double _speed_howManyAdsorption;
	static double _speed_howManyDesorption;
	static double _speed_howManyPrecursor;
	static float    _speed_rateFactor;
	static float  _speed_oldkt;
};

} /* namespace OKMC */
#endif /* LATTICEATOMDEPI_H_ */
