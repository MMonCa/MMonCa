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
//-------------------------------------------------------------------------
// DEFINES
//-------------------------------------------------------------------------
#ifndef _INCLUDE_DEFS_H_
#define _INCLUDE_DEFS_H_

#ifndef UNIX
 // No escapes ANSI
 #define _CLS_    ""
 #define _NORMAL_ ""
 #define _BOLD_   ""
 #define _RED_    ""
 #define _GREEN_  ""
 #define _BROWN_  ""
 #define _BLUE_   ""
 #define _OK_     ""
 #define _ERROR_  " --> ERROR"

#else
#ifndef NO_ANSI
 #define _CLS_   "\033[H\033[2J"
 #define _NORMAL_ "\033[0m"
 #define _BOLD_   "\033[1m"
 #define _RED_    "\033[31m"
 #define _GREEN_  "\033[32m"
 #define _BROWN_  "\033[33m"
 #define _BLUE_   "\033[34m"
 #define _OK_     "\033[70G\033[32m\033[1m CORRECT\033[0m"
 #define _ERROR_  "\033[70G\033[31m\033[1m ERROR\033[0m"
#else
 #define _CLS_    ""
 #define _NORMAL_ ""
 #define _BOLD_   ""
 #define _RED_    ""
 #define _GREEN_  ""
 #define _BROWN_  ""
 #define _BLUE_   ""
 #define _OK_     ""
 #define _ERROR_  " --> ERROR"
#endif

#endif

#define READ_ERROR                  300
#define WRITE_ERROR                 301
#define MEMORY_ERROR                302

#define NOATOMS_ERROR               400
#define MAX_LATTICE_SITES_ERROR     401
#define RARE_EVENT_ERROR            402
#define HISTOGRAM_ERROR             403
#define MAX_ATOMS_ERROR             404
#define NO_LAYERS_ERROR             405
#define NO_SIM_ERROR                406
#define OUT_OF_RANGE_ERROR          407

#define FLAT_NOT_DEFINED_ERROR      500
#define DOSE_DISTRIBUTION_ERROR     501
#define ELECTRONIC_STOPPING_ERROR   502
#define IONIZATION_ERROR            503

#define INCOMPLETE_STRING_ERROR     600
#define NOT_END_WITH_ZERO_ERROR     601
#define BOOLEAN_ERROR               602

// Used for comparisons of values in many places in the program
//  In lattice units (1.0862e-5 )
#define DELTA1  1.0862e-5

// Define number of entries in Therm.Table
#define NTBL 200

// Define the maximum number of lattice sites
#define NEIGHMAX 30

// Define the maximum number of atoms involved in the simulation
#define N_ATOMS 10

// Defines for List<>.Sort function
#define sortASCENDENT   1
#define sortDESCENDENT -1

// Conversion factor from Degrees to Radians.
#define DEGRAD   0.017453292519943296

#define TRANSLATION_PARAMETER   0.5

#define NX 112
#define NY 112
#define NZ 112

#define BODY_MAX 20

// Using the M.K.S. units system
#define proton_mass     1.672623e-27   /* Proton mass     (Kg)  */
#define electron_charge 1.60217733e-19 /* Electron charge (C)   */
#define bohr_radius     5.29177249e-11 /* Bohr radius     (m)   */
#define bohr_velocity   2.18769e+06    /* Bohr velocity   (m/s) */

#define CTE 2*electron_charge/(proton_mass*bohr_velocity*bohr_velocity)

#ifdef INTEL_C
#define M_PI 3.141592654
#define M_E 2.718281828
#endif

// Number of elements of the probability table for thermal amplitudes
#define THERMAL_NTBL  2048 // 2^{11}

#define TABLEDIRDEFAULT "./"
#define EDTDIRDEFAULT   "./"

#endif //_INCLUDE_DEFS_H_
