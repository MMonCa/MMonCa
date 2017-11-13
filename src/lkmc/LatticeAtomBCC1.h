/*
 * LatticeAtomBCC1.h
 *
 *  Created on: Dec 19, 2013
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


#ifndef LATTICEATOMBCC1_H_
#define LATTICEATOMBCC1_H_

#include "LatticeAtom.h"

namespace LKMC {

class LatticeAtomBCC1 : public LatticeAtom
{
public:
	LatticeAtomBCC1(LASTATE, Kernel::Domain *p, Kernel::M_TYPE, Kernel::P_TYPE, const Kernel::Coordinates &c);
	LatticeAtomBCC1(std::istream &is);
	virtual ~LatticeAtomBCC1() {}

	virtual LatticeParam::TYPE getLatticeType() const { return LatticeParam::BCC; }
	virtual void perform (Kernel::SubDomain *, unsigned eventType);
	virtual float getRate(unsigned eventType, float kT) const;
	virtual void setIndex(unsigned eventType, int idx) { _idx = idx; }
	virtual int  getIndex(unsigned eventType) const { return _idx; }
	virtual unsigned getNEvents() const { return 1; } //recrystallization

	//              4            22           70         Perfect, but due to periodic BC it can be more.
	enum { FIRSTN = 8, SECONDN = 16, THIRDN = 24 };
	virtual unsigned first() const  { return FIRSTN; }
	virtual unsigned second() const { return SECONDN; }
	virtual unsigned third() const  { return THIRDN; }
	virtual void  setNeiDist(int i, float d2) {  }
	virtual float getNeiDist(int i) const { return 0; }
	virtual void dumpXYZ(const std::vector<std::pair<int, Kernel::Coordinates> >&, const std::string &name) const;
	virtual bool isDefective() const { return false; }
	virtual std::string getClass() const { return "BCC1"; }

	virtual void restart(std::ostream &) const;

private:
	int _idx;
};

} /* namespace LKMC */
#endif /* LATTICEATOMBCC1_H_ */
