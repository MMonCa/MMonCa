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

#ifndef KMCMCCLIENT_H
#define KMCMCCLIENT_H

#include "kernel/Coordinates.h"
#include "kernel/Material.h"
#include "kernel/ParticleType.h"
#include "io/GetMaterial.h"
#include "io/MeshParser.h"
#include <map>
#include <string>
#include <tcl.h>
#include "kernel/RNG.h"

namespace IO { class ParameterManager; }

namespace OKMC { class Defect; class Particle; class Cluster; class Interface; }
namespace Kernel { class MeshElement; class SubDomain; class Domain; }

namespace Domains {

class MCClient
{
public:
    MCClient(Tcl_Interp *p, const Kernel::Coordinates &m, const Kernel::Coordinates &M,
    	const std::string &proc);
    MCClient(std::istream &);
    MCClient(const IO::MeshParser *);

    ~MCClient();
    
    void beginInsert();
    void endInsert();
    OKMC::Cluster *   createMC(const std::string &name, const std::string &type, const Kernel::Coordinates &c);
    OKMC::Defect *    createMP(const std::string &name, const std::string &stat, const Kernel::Coordinates &c, bool bReact);
    OKMC::Interface * createIF(const std::string &name, const Kernel::Coordinates &c);

    const IO::GetMaterial * getMaterial() const { return _getMaterial; }
    void resetGetMaterial();  //to be called when the information from getMaterial is no longer available

    void getCell(Kernel::Coordinates &m, Kernel::Coordinates &M) const { m=_min; M=_max; }
    bool isInCell(const Kernel::Coordinates &c) const {return c.isInto(_min, _max); }
    double rand() { return _rng2.rand(); }

    void checkAmorphization(Kernel::SubDomain* pSub, const Kernel::Coordinates &c);
    void checkLocalAmorphization(Kernel::MeshElement *,std::vector<const Kernel::MeshElement *> );
	std::map<std::string, unsigned> getCreatedDefects(Kernel::M_TYPE mt) const {return _created[mt];}

    void restart(std::ostream &) const;
    bool isFromStream() const { return _isFromStream; }

private:
    Kernel::Coordinates _min, _max;
    Kernel::RNG _rng2; //for Interfaces
	IO::GetMaterial *_getMaterial;

	std::map<std::string, unsigned> _created[Kernel::MAX_MATERIALS];
	std::map<std::string, unsigned> _discarded;
	std::map<std::string, unsigned> _out;

	const bool _isFromStream;
};

}

#endif
