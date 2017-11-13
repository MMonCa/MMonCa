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

#ifndef KMCEVENT_H
#define KMCEVENT_H

#include <iostream>

namespace Kernel
{
class Domain;
class SubDomain;

class Event
{
public:
    enum E_TYPE { LATTICEATOM, MOBILEPARTICLE, INTERFACE, CLUSTER, EMPTY, SINK, CASCADE, CELL, UNDEFINED_EVENT };

    Event(Domain *p) : _pDomain(p) {}
    Event(std::istream &);
    
    virtual             ~Event() {}
    virtual void        perform (Kernel::SubDomain *, unsigned eventType) = 0;
    virtual void        setIndex(unsigned eventType, int idx) = 0;
    virtual unsigned    getNEvents() const = 0;
    virtual int         getIndex(unsigned eventType) const = 0;
    virtual float       getRate(unsigned eventType, float kT) const = 0;
    virtual E_TYPE      getEType() const = 0;

    static const char * getEName(E_TYPE ev)  { return _names[ev]; }
    const Domain *      getDomain() const    { return _pDomain; }

    void restart(std::ostream &) const;

protected:
    Domain *_pDomain;
private:
    static const char *_names[UNDEFINED_EVENT];
};

}

#endif
