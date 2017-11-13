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
// MAIN.H
//
///////////////////////////////////////////////////////////////////////////
//#define NDEBUG
 
//-------------------------------------------------------------------------
// INCLUDES
//-------------------------------------------------------------------------
#include "defs.h"
#include "vector.h"
#include "list.h"
#include "mainloop.h"

#ifdef DAMAGE1D
 #include "amorf.h"
#endif
#ifdef DAMAGE3D
 #include "amorf3d.h"
#endif

namespace BCA {
//-------------------------------------------------------------------------
// GLOBAL VARIABLES
//-------------------------------------------------------------------------
char DIRECTORY[300];
short unsigned int SpecificPotential[100][100];
int LayerNumber;
double RI = 3.0; // 3 Angstroms
List<Vector> LP;

List<AtomDefinition>       Atoms;    // Atoms used in simulation
List<LayerDefinition>      Bulk;
List<SimulationDefinition> SIMULACIONES;

long double NSimMean;
long double ApsisIterations;   // For statistics on apsis calculus.
long PEPE=0;        // Counter for simultaneous collisions done
int PROJS[N_ATOMS];
}
