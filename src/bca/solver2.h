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
// SOLVER2.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_SOLVER2_H_
#define _INCLUDE_SOLVER2_H_

#include "table3d.h"
namespace BCA {
///////////////////////////////////////////////////////////////////////////
// SOLVER class
class Solver
{
 private:
  // Variables ------------------------------------------------------------ 
      double EP, Ghi, p;
    Vector RT, RP, VP, DeltaR;
       int AtomP, AtomT;
       
 public:
      double Q;
       TBL Table; 

  // Functions ------------------------------------------------------------
    void FirstInit( );
    void All ( Projectile &P, Projectile &T, double &AdvanceMinimum,
          double &Einel, double &QElect );
};
}
#endif  // _INCLUDE_SOLVER2_H_
