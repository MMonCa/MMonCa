/*
 * LatticeAtomDSPER.h
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

#ifndef LATTICEATOMDSPER_H_
#define LATTICEATOMDSPER_H_

#include "LatticeAtomDiamond.h"

namespace LKMC {

class LatticeAtomDSPER: public LKMC::LatticeAtomDiamond
{
public:
	LatticeAtomDSPER(LASTATE, Kernel::Domain *p, Kernel::M_TYPE, Kernel::P_TYPE, const Kernel::Coordinates &c, unsigned orientation);
	LatticeAtomDSPER(std::istream &);
	virtual ~LatticeAtomDSPER() {}

	virtual void perform (Kernel::SubDomain *, unsigned eventType);
	virtual float getRate(unsigned eventType, float kT) const;

	virtual void setIndex(unsigned eventType, int idx) { _idx = idx; }
	virtual int  getIndex(unsigned eventType) const { return _idx; }
	virtual unsigned getNEvents() const { return 1; } //recrystallization

	//              4         4+12=16   16+12=28         Perfect, but due to periodic BC and twins it can be more.
	enum { FIRSTN = 9, SECONDN = 30, THIRDN = 64 };
	virtual unsigned first() const  { return FIRSTN; }
	virtual unsigned second() const { return SECONDN; }
	virtual unsigned third() const  { return THIRDN; }

	virtual std::string getClass() const { return "DSPE"; }
	virtual void restart(std::ostream &) const;

private:
	int _idx;

	void event100(Kernel::SubDomain *, std::map<int, LatticeAtomDiamond *> &);
	void event110(Kernel::SubDomain *, std::map<int, LatticeAtomDiamond *> &);
	void event111(Kernel::SubDomain *, std::map<int, LatticeAtomDiamond *> &);
	float getSPERRate(float kT) const;
	unsigned getPlane(unsigned coord) const;

	void     direction100(ublas::vector<double> &) const;
	void     direction110(ublas::vector<double> &) const;
	double   stressCorrection_2ndOrder(unsigned, unsigned, float) const;

	void     cross_product(const ublas::vector<double> &, const ublas::vector<double> &, ublas::vector<double> &) const;
	unsigned idxMinVector(std::vector<double> &) const;
	unsigned idxMaxVector(std::vector<double> &) const;
	void     normalize(ublas::vector<double> &) const;

	double getGFLSFactor();

	mutable int _plane;
};

} /* namespace LKMC */
#endif /* LATTICEATOMDSPER_H_ */
