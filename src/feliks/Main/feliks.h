/* Finite Element Method Module
 *
 * Author: ignacio.romero@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain,
 *      and       Technical University of Madrid (UPM), Madrid, Spain
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
/*
 *  feliks.h
 *
 *  F EL I K S,    a  general purpose finite element code
 *  by Ignacio Romero
 *  Universidad Polit�cnica de Madrid.
 *
 *  Created by ignacio on Sun Mar 02 2003.
 *  Last modification, Dec 21 2010
 *
 *
 *  General options and defintions for FELIKS.
 *  This options are valid only at compilation time: after modifications, the
 *  program needs to be recompiled
 *
 */

#ifndef _feliks_h
#define _feliks_h

#define  FELIKSPATH       "/Users/ignacio/Programs/feliks/Source"

// FELIKS built information
#define  FELIKS_VERSION          3.40
#define  FELIKS_DATE             "August 2012"


// these constants define the default convergence criteria in Newton-Raphson iterations
// they can be modified in the input file 
#define  RESIDUAL_TOLERANCE    1e-9
#define  ENERGY_TOLERANCE      1e-16    // Relative error. This is the main convergence tolerance
#define  ENERGY_ABS_TOLERANCE  1e-20
#define  MAX_ITERATIONS        10       // max. # of iterations in a NR loop
#define  MAX_NRTRIALS          4        // max. # of NR loops before declaring an error


// these variables set the maximum sizes of elements and help to speed up memory allocations
#define  FELIKS_MAX_ELNODES      10
#define  FELIKS_MAX_FACENODES		9
#define  FELIKS_MAX_SURFACEGP		9

#define TETLIBRARY

#endif

/*! \mainpage F EL I K S, a general purpose finite element program for research...
 *
 */
