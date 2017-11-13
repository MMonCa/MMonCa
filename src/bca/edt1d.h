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
///////////////////////////////////////////////////////////////////////////
//
// EDT1D.H
//
// Version 2.2a October 2001
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_EDT1D_H_
#define _INCLUDE_EDT1D_H_

// #define NOEXTRAPOLATION
//  If not defined makes a extrapolation else no extrapolation
//   in Interpolate1D() and GetRs()

#include "defs.h"
#include "indexes.h"

#include <stdlib.h>
#include <fstream>
#include <string.h>
#include <math.h>

#define NMXV 2000
namespace BCA {
// Class: Electronic Density Table 1D
class EDT1D
{
  protected:
    // EDT simple
    double Radius [ NMXV ];
    double Density[ NMXV ];
    int    ND;   // Number of Data

  public:
    EDT1D();
    void   ReadTable(const std::string &FileName, int PP=1 );
    double GetMaxR();
    double Interpolate1D( double R );
    double GetRs( double r );
};
}
#endif	// _INCLUDE_EDT1D_H_
