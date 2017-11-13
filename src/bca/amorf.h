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
//  AMORF.H
//  Modified Kinchin-Pease model for Damage Accumulation
//
////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_AMORF_H_
#define _INCLUDE_AMORF_H_

#include "defs.h"
#include "varhst.h"
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <string>
#include <stdlib.h>

namespace BCA {

extern char DIRECTORY[];

class clase_Amorfizacion
{
 private:
#ifdef FAST_KP 
  double KK3;
#else
  double KK, KK2;
#endif  
 public:
  VarHst DefectDensity;     // The positions are in (A)
#ifdef FAST_KP
  VarHst E_Acc;
#else  
  VarHst AA, BB;
#endif  
  double Area;              // (A^2)
  double CutOffEnergy;      // Cut-off energy by volume unit (eV/A^3)
  double DisplacementEnergy;// Displacement energy (eV)
  double frec;              // Surviving recombination factor
  double k;                 // Kinchin-Pease constant
  double Na;                // Amorphization density (defects/A^3)

  int LevelOfSimulation;
  int NNN;

  double DoseRate;
#ifdef MODEL2
  double CteA;
  double CteB;
  double CteC;
  double CteD;
  double CteK;
  double ancho;
  double Weff;
  
  double Mean_CrossSection;
  double TotalGeneratedFrenkelPairs;

  double DeltaDamageByCascade;
#endif
#ifdef ENERGIA_DEPOSITADA 
  double ETOTAL;
  ofstream OO;
#endif 

  //  Constructor ---------------------------------------------------------
  clase_Amorfizacion();

  // Init -----------------------------------------------------------------
  void Init( double Dose, double IA, int LS,
        double RecombinationFactor, double KinchinPeaseConstant,
        double Na_, double COE, double DE=15.0,
        int MultipleImplant=0, const std::string &PreviousDamageFile="",
    double DR=1e13, double WE=1.0 );

  // Damage accumulation --------------------------------------------------
  void DamageAccumulation
       ( double posx, double energia, char tipo,double Weigth );

  // Update defect density ------------------------------------------------
  void UpdateDefectDensity();

  // Damage accumulation. Sequential version ------------------------------
  void SeqDamageAccumulation
       ( double posx, double energia, char tipo='I',double Weigth=1.0, double Temp=300 );

  // Write amorphization histogram ----------------------------------------
  void WriteAmorphization();

  // SUMA Energia ----------------------------------------------------------
  #ifdef ENERGIA_DEPOSITADA
  void Suma(double X, double Y, double Z, double Energy,
               char C[2], double Weigth=1);
  #endif

  // Destructor ------------------------------------------------------------
  ~clase_Amorfizacion();

};

void Plot2PSDamage( FILE * OutPS, const std::string &File );

}

#endif //_INCLUDE_AMORF_H_
// ------------------------------------------------------------------------
