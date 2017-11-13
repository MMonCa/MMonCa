 /*
 * Original author: jesus.hernandez.mangas@tel.uva.es
 *
 * Copyright 2014 University of Valladolid, Valladolid, Spain.
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
 *
 * Adapted to MMonCa: I. Martin-Bragado. IMDEA Materials Institute.
 *
 */ 
//////////////////////////////////////////////////////////////////////////
//
// THERMAL.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_THERMAL_H_
#define _INCLUDE_THERMAL_H_

#include "defs.h"
#include "vector.h"
#include "list.h"
#include "atomdef.h"
namespace BCA {
class ThermClass
{
  public:
    double Table[THERMAL_NTBL];   // Probability table
    double TAmp[103]; // Thermal amplitudes

    void Generate( List<AtomDefinition> &Atoms, double T );
};
}
#endif  // _INCLUDE_THERMAL_H_
