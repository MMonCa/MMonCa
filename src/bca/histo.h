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
// HISTO.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_HISTO_H_
#define _INCLUDE_HISTO_H_

#include <stdio.h>
#include <fstream>
#include <string.h>

#include "vector.h"
#include "list.h"
#include "defs.h"
#include "mainloop.h"

#define LOG2_10 3.321928095
#define N_BIN 100
#define N_BIN2D 50
#define MAX_SPLIT 20
namespace BCA {
extern char DIRECTORY[300];
// ------------------------------------------------------------------------
class VirtualIon
{
 public:
  double Depth;
  double Weigth;
  double Other;
};

class Hist
{
 private:
  double Dose;
  double RealIons;
  double DMaxO;
  List<VirtualIon> LVI; // Stopped virtual ions

 public:
  double H[N_BIN+1];      // Histogram for depths
  double HR[N_BIN+1];     // Histogram for distances 
  double H_err[N_BIN+1];  // Absolute error
  double DMax, DMin;
  double Outside;         // Backscatered/sputtered ions
  
  void Initialize(double MaxDepth=1000.0/*(A)*/,double D=1e13 );
  void Add( double D, double W, double O );
  void Do();         
  void SaveHisto (int NProj, int HST,char* HSTFile,
                    int FullDoseOutput, int KK, unsigned long TotalRealIons );
  void SaveHistoR(int NProj, int HST,char* HSTFile);
  void ShowFinalInfo(int NProj);
  void SolvePearsonIV(int NProj,char * HSTFile );
};

// ========================================================================
// New class to describe a complete virtual ion ///////////////////////////
class IonVirtual
{
 public:
    Vector R;
    double Weigth;
    double Other;
};

// New class to solve a 3D histogram //////////////////////////////////////
class Histo3D {
 private:
  List<IonVirtual> LIV;
 public:
  double H3D[N_BIN2D][N_BIN2D][N_BIN2D]; // Locates at the stack

  void Initialize(){ LIV.DeleteAll(); }
  void Add( Vector R, double W, double O );

  void Solve3D( int NProj, double Dose, double WinAB, int ascii=0 );
  void Show( int NProj, double Dose );
};
}
#endif  // _INCLUDE_HISTO_H_
