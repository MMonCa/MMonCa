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
// SIMUL.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_SIMUL_H_
#define _INCLUDE_SIMUL_H_

#include "iostream"
#include "vector.h"
#include "stopping.h"
namespace BCA {
class SimulationDefinition
{
 public:
    // Variables
        int AtomP, PROJS[N_ATOMS], Miller;
        double ENERGY, Tha, Phi, Divergency, Dose, Temperature, RS0[MAXLAYER];
        Vector DIR;
        unsigned long NumberOfImplants;
        double DoseDistribution[10];
        int  DoseSplitting;
    // CONSTRUCTOR
    SimulationDefinition()
    {
        Default();
    };
    // Default values
    void Default()
    {
        AtomP=1;
        for(int i_=0; i_<N_ATOMS; i_++) PROJS[i_]=0;
        ENERGY=200.0;
        Miller=0;
        DIR(1,0,0);
        Tha=0.0;
        Phi=0.0;
        Divergency=0.0;
        NumberOfImplants=10;
        Dose=1.0e13;
        Temperature=300.0;
        for( int k=0;k<MAXLAYER;k++ ) RS0[k]=1.85;
        DoseSplitting = 0;
        for( int j=0; j<10; j++ ) DoseDistribution[j]=0.0;
    };
    // Output stream
    friend std::ostream &operator <<( std::ostream &Out, SimulationDefinition A );
};
}
#endif  // _INCLUDE_SIMUL_H_
