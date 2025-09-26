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

#ifndef LKMCLATTICEATOM_H
#define LKMCLATTICEATOM_H

#include "kernel/Coordinates.h"
#include "kernel/Event.h"
#include "LatticeSite.h"
#include "Lattice.h"
#include <vector>

namespace Kernel
{
	class MeshElement;
	class Domain;
}

namespace LKMC
{

class LatticeAtom : public Kernel::Event, public LatticeSite
{
private:
    static constexpr int NEIGH_TOO_CLOSE = -1;
    static constexpr int NEIGH_NO_PLACE  = -2;
    static constexpr float MIN_DIST_FACTOR = 0.8;

public:
	enum LASTATE { LS_AVAILABLE, LS_PERFORMED, LS_PRECURSOR };

	LatticeAtom(LASTATE, Kernel::Domain *p, Kernel::M_TYPE, Kernel::P_TYPE type, const Kernel::Coordinates &c);
	LatticeAtom(std::istream &);
    virtual ~LatticeAtom();

    bool getPerformed() const { return _state == LS_PERFORMED; }
    void displaceCoordinates(Kernel::Coordinates &);
    void setState(LASTATE st);
    LASTATE getState() const { return _state; }
    
    virtual std::string getClass() const = 0;
    virtual LatticeParam::TYPE getLatticeType() const = 0;
    virtual void perform (Kernel::SubDomain *, unsigned eventType) = 0;
    virtual void setIndex(unsigned eventType, int idx) = 0;
    virtual int  getIndex(unsigned eventType) const = 0;
    virtual unsigned getNEvents() const = 0;
    virtual float getRate(unsigned eventType, float kT) const = 0;
    virtual Event::E_TYPE getEType() const { return LATTICEATOM; }
    virtual bool isDefective() const = 0;

    LatticeAtom * getNeighbor(unsigned i) const { return _neighbors[i]; }
    int getNumber() const { return _number; }
    int getCoordination(unsigned i) const;

    virtual unsigned first()  const = 0;
    virtual unsigned second() const = 0;
    virtual unsigned third()  const = 0;
    virtual void  setNeiDist(int, float) = 0;
    virtual float getNeiDist(int) const = 0;

    virtual void restart(std::ostream &) const;

protected:
	static int _num;
	int _number;
	LASTATE _state;
	LatticeAtom **_neighbors;
	unsigned char _maxNeigh;
	
    int tryNeig(unsigned, LatticeAtom *pLA, float dist2, float minDist2);
    void insertNeig(LatticeAtom * const pLA, int const where, float dist2);
	void insertNeighbors();
	void removeNeighbor(LatticeAtom *);
	void updateME4Epitaxy(Kernel::SubDomain *); //to add more "non-crystalline" atoms...
};

}

#endif
