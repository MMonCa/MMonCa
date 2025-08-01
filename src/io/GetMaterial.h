/*
 * GetMaterial.h
 *
 * Created on: Jun 20, 2012
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

#ifndef GETMATERIAL_H_
#define GETMATERIAL_H_

#include "kernel/Coordinates.h"
#include "kernel/Material.h"
#include <cstdint>
#include <vector>
#include <tcl.h>

namespace IO {

class GetMaterial {
public:
	GetMaterial(Tcl_Interp *p, const std::string &procName) : _pTcl(p), _procName(procName) {}
	GetMaterial(std::istream &is);
	
	void fill(uint32_t const aDomain, uint32_t const aLinesXsize, uint32_t const aLinesYsize, uint32_t const aLinesZsize,
        uint32_t const aLinesZindexMin, uint32_t const aLinesZindexMax, std::vector<Kernel::M_TYPE> const * const aMaterials);

	Kernel::M_TYPE operator()(const Kernel::Coordinates &c, unsigned idx) const;
	static void restart(std::ostream &);

private:
	Tcl_Interp *_pTcl;
	std::string _procName;  //when using a procedure
	std::vector<std::vector<Kernel::M_TYPE> >_materials; //when reading from file

};

} /* namespace IO */
#endif /* GETMATERIAL_H_ */
