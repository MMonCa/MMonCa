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
////////////////////////////////////////////////////////////////////////////
//
// SSF.H    Specific Screening Function
//
////////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_SSF_H_
#define _INCLUDE_SSF_H_

#include "defs.h"
#include "vector.h"

#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <string.h>

#define N_SSF 1000
#define D_SSF 5.0

#define ZBL_SCREENING   0
#define SPECIFIC        1
#define THOMAS_FERMI    2
#define MOLIERE         3
#define LENZ_JENSEN     4
#define BOHR            5
#define SPECIFIC_WOPD   6
namespace BCA {
extern short unsigned int SpecificPotential[100][100];

////////////////////////////////////////////////////////////////////////////

void GetSSFName(int nZ1, int nZ2, char *Name );

////////////////////////////////////////////////////////////////////////////

class SSF
{
 public:
	double _au;

    void Init( int Z1, int Z2, int Tipo );

    void ZBL_Screening( double X, double &S, double &DS );
    void Thomas_Fermi_Screening( double X, double &S, double &DS );    
    void Moliere_Screening(double X, double &S, double &DS );    
    void Lenz_Jensen_Screening(double X, double &S, double &DS );
    void Bohr_Screening(double X, double &S, double &DS );

    double Screening( double R, double &S, double &DS );    
    double ScreeningOnly( double R, double &S );
 private:
    double _screen[N_SSF];
    double _dscreen[N_SSF];

};

}
#endif  // _INCLUDE_SSF_H_
