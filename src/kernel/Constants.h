/*
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain
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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#define PLANCK         6.626176e-34      // J.s
#define PLANCK_BAR     1.054571726e-34   // J.s
#define M0             9.10938215e-31    // kg
#define KB             1.3806488e-23     // J/K
#define Q              1.60217657e-19    // C
#define AVOGADRO       6.02214129e23     // /mol
#define BOHR_RADIUS    5.2917721092e-11  // m
#define COULOMB        8.9875517873681e9 // N m^2 C^-2

#define CELSIUS_TO_KELVIN(x)        ((x) + 273.15)
#define KELVIN_TO_CELSIUS(x)        ((x) - 273.15)

#define ELECTRONVOLT_TO_HARTREE(x)  ((x) / 27.21138386)
#define HARTREE_TO_ELECTRONVOLT(x)  ((x) * 27.21138386)

#endif /* ! CONSTANTS_H */
