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

#ifndef KMCRATEMANAGER_H
#define KMCRATEMANAGER_H

#include "Event.h"
#include <fstream>
#include <string>
#include <vector>

namespace LKMC { class LatticeAtom; }
namespace MMC  { class MMCAtom;     }
struct Tcl_Interp;
namespace Kernel
{
class Domain;
class SubDomain;
class MeshElement;

class RateManager
{
public:
    RateManager(Domain *);
    ~RateManager();

    void setTempK(float K);
    void anneal(double endTime, bool bDepth, float depth, long unsigned events);
    void setDepthLA(float depth);
    float getDepthLA() const { return _depthLA; }
    long unsigned getEvents() const { return _nEvents; }
    double getTime() const { return _time; }
    void   setTime(double t) { _time = t; }
    float  getkT() const { return _kT; }
    float  getT() const { return _kelvin; }

    //MeshElement needed to compute the subdomain and level
    void insert(Event *, MeshElement *);
   	void update(Event *, MeshElement *);
    void remove(Event *, MeshElement *);

    void printOutput(double printtime) const;

    SubDomain * getSubDomain(unsigned i) const { return _subDomains[i]; }
    unsigned    getSubDomains() const { return _subDomains.size(); }

    void restart(std::ostream &) const;
    void restart(std::istream &);

private:
    Domain *_pDomain;
    std::vector<SubDomain *> _subDomains;
    float _kB;
    float _kelvin;
    float _kT;
    long unsigned _nEvents;
    float _depthLA;
    double _time, _timeLastEvent, _timeNextEvent;
    double _lastMaxRate;
    bool _bAveraged;
    unsigned _nLevels;
};

}

#endif
