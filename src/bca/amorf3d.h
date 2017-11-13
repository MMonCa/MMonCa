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
////////////////////////////////////////////////////////////////////////
//
//  AMORF3D.H
//  Modified Kinchin-Pease model for Damage Accumulation
//
////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_AMORF_H_
#define _INCLUDE_AMORF_H_

#include "varhst.h"
#include "defs.h"
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <stdlib.h>

namespace BCA {

extern char DIRECTORY[];

class clase_Amorfizacion3D
{
 private:
#ifdef FAST_KP 
  double KK3;
#else
  double KK, KK2;
#endif  
 public:
  VarHst3D DefectDensity;     // The positions are in (A)
#ifdef FAST_KP
  VarHst3D E_Acc;
#else  
  VarHst3D AA, BB;
#endif  
  double Area;              // (A^2)
  double WinA, WinB;        // (A)
  double min_y,max_y,min_z,max_z; // Indexes
  double CutOffEnergy;      // Cut-off energy by volume unit (eV/A^3)
  double DisplacementEnergy;// Displacement energy (eV)
  double frec;              // Surviving recombination factor
  double k;                 // Kinchin-Pease constant
  double Na;                // Amorphization density (defects/A^3)

  int LevelOfSimulation;
  int NNN;

#ifdef MODEL2
  double DeltaDamageByCascade;
#endif
#ifdef ENERGIA_DEPOSITADA
  double ETOTAL;
  ofstream OO; 
#endif 

  //  Constructor ---------------------------------------------------------
  clase_Amorfizacion3D();

  // Init -----------------------------------------------------------------
  void Init3D( double Dose, double IA, int LS,
        double RecombinationFactor, double KinchinPeaseConstant,
        double Na_, double COE, double WA, double WB, double DE=15.0,
        int MultipleImplant=0, const std::string &PreviousDamageFile="" );

  // Damage accumulation --------------------------------------------------
  void DamageAccumulation3D
       ( Vector pos, double energia, char tipo, double Weigth );

  // Update defect density ------------------------------------------------
  void UpdateDefectDensity3D(Vector InitPos);

  // Damage accumulation. Sequential version ------------------------------
  void SeqDamageAccumulation3D
       ( Vector pos, double energia, 
            char tipo='I',double Weigth=1.0, double Temp=300.0 );

  // Write amorphization histogram ----------------------------------------
  void WriteAmorphization3D(int Opciones);

  // Add Energy (test)----------------------------------------------------------
  #ifdef ENERGIA_DEPOSITADA
  void Suma(double X, double Y, double Z, double Energy,
               char C[2], double Weigth=1);
  #endif

  // Destructor ------------------------------------------------------------
  ~clase_Amorfizacion3D();

};
}
#endif //_INCLUDE_AMORF_H_
 // ------------------------------------------------------------------------
