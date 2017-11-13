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
// STOPPING.H   
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_STOPPING_H_
#define _INCLUDE_STOPPING_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <string.h>
#include "defs.h"

#define NTABLEX 300
#define NTABLEY 2750
#define DLIMX   6.0
#define DLIMY   5.5

#define MAXLAYER 5
namespace BCA {
double ZBLStopping ( double v, double z, double rs );
double OurStopping( double v, double z, double rs );
double NoStopping  ( double v, double z, double rs );

double ZBL_I( double yr ); 
double  BK_I( double yr ); 
double CGJ_I( double yr ); 
double  MP_I( double yr ); 

struct StoppingStruct
 {
  char Name[10];
  char Comment[80];
  double (*Stopping)( double v, double z, double rs );
 };
 
struct IonizationStruct 
 {
  char Name[10];
  char Comment[80];
  double (*Ionization)( double yr );
 };

class STP
{
 public:
    double * TABLA[104][MAXLAYER]; 
 
    STP();
    void Clear();
    double Stopping( double V, double dZ, double rS );
};
}
#endif  // _INCLUDE_STOPPING_H_
