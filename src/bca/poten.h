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
// POTEN.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_POTEN_H_
#define _INCLUDE_POTEN_H_

#include "vector.h"
#define NUMMAX 15  // Max value cached
#define MATR 2000  // Number of sampled points to cache
namespace BCA {

class ZBL
{
 private:
    #ifdef CACHED_ZBL
    double TABLEV [MATR];     // Cache
    double TABLEDV[MATR];
    double DX;
    #endif
 public:
    double EXPAR[ 4 ];        // Coefs.
    double  VPAR[ 4 ];
    double DVPAR[ 4 ];

    ZBL();
    void Potential( double Xa, double &VR, double &DVR );
    void PotentialOnly( double Xa, double &VR );
};
}
#endif  // _INCLUDE_POTEN_H_
