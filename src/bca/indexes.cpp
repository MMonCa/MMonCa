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
// INDEXES.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "indexes.h"
namespace BCA {
double E_GRID[N_ENERGY] = { 
    10.,    12.11527659, 14.67799268, 20.08546611, 21.5443469, 26.10157216, 
    31.6227766, 38.31186849, 46.41588834, 56.23413252, 68.12920691, 82.54041852, 
   100.,   121.1527569, 146.7799268, 200.8546611, 215.443469, 261.0157216, 
   316.227766, 383.1186849, 464.1588834, 562.3413252, 681.2920691, 825.4041852,
  1000.,  1211.527569, 1467.799268, 2008.546611, 2154.43469, 2610.157216, 
  3162.27766, 3831.186849, 4641.588834, 5623.413252, 6812.920691, 8254.041852,
 10000., 12115.27569, 14677.99268, 20085.46611, 21544.3469, 26101.57216, 
 31622.7766, 38311.86849, 46415.88834, 56234.13252, 68129.20691, 82540.41852,
100000.,121152.7569, 146779.9268, 200854.6611, 215443.469, 261015.7216, 
316227.766, 383118.6849, 464158.8834, 562341.3252, 681292.0691, 825404.1852,
1000000.,1211527.569, 1467799.268, 2008546.611, 2154434.69, 2610157.216, 
3162277.66, 3831186.849, 4641588.834, 5623413.252, 6812920.691, 8254041.852
 };

double S_GRID[N_S] = { 0.0005, 0.001, 0.005, 0.01, 0.03, 0.05, 
                       0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 
                       0.8, 1, 1.2, 1.4, 1.6, 1.8, 
                       2, 2.2, 2.4, 2.6, 2.8, 3, 
                       3.2, 3.4, 3.6, 3.8, 4.0, 5.0 };
}
