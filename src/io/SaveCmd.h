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

#ifndef GENERALSAVECMD_H
#define GENERALSAVECMD_H

#include "Command.h"
#include "domains/SimData.h"

namespace Kernel
{
	class Mesh;
}
namespace IO {


class SaveCmd : public Command
{
public:
    SaveCmd(Tcl_Interp *, int argc, const char *argv[]);
    ~SaveCmd() {};
    virtual int operator()();

    enum FORMAT_TYPES { XYZ_TYPE, LAMMPS_TYPE, ATOMEYE_TYPE, CSV_TYPE, VTK_TYPE };
    void save(const std::string &filename, FORMAT_TYPES, bool bDefects, const std::vector<std::string> &defects) const;

private:
    void writeAtomeye(std::ofstream &, std::vector<Domains::ParticleData> &) const;
};

}

#endif
