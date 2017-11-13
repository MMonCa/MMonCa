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
//////////////////////////////////////////////////////////////////////////
//
//  TABLE3D.H 
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_TABLE3D_H_
#define _INCLUDE_TABLE3D_H_

#define HEADERSIZE      200

#include "atom.h"
#include "indexes.h"
#include "col3d.h"
namespace BCA {
extern List<AtomDefinition> Atoms;

///////////////////////////////////////////////////////////////////////////
enum TableVariable
 {
  _X1 = 0,  /* X1   Value (A)  */   // Output variables
  _X2,      /* X2   Value (A)  */
  _COS2,    /* Cos2 Value ()   */
  _E1,      /* E1   Value (eV) */
  _E2,      /* E2   Value (eV) */
  _Q,       /* Q    Value (eV) */
  _QElect   /* QElect Value (eV) */
 };

typedef struct  {
    double TABLE[ N_ENERGY ][ N_PHI ][ N_THETA ][ N_ALPHA ][ N_S ][ NV3D ];
} tTABLE;

///////////////////////////////////////////////////////////////////////////

// Define the TaBLe class /////////////////////////////////////////////////
class TBL
{
 private:
    tTABLE * GRAN[ N_ATOMS ][ N_ATOMS ];
    double rs0_table;
    char* TABLEDIR;
    char Cabecera[HEADERSIZE+1];
 public:
    // Constructor
    TBL();
//  void GetName( int nZ1, int nZ2, char *Name );

    // Read Tables from Disk
    void ReadTable ();
    int  LoadASCII ( int nZ1, int nZ2, char * FileName );
    int  LoadBINARY( int nZ1, int nZ2, char * FileName );
    void SaveASCII ( int nZ1, int nZ2, char * FileName );
    void SaveBINARY( int nZ1, int nZ2, char * FileName );

    void CalculeTable( int A1, int A2 );

    // Interpolate variables
    void InterpolateAll_( int Atom1, int Atom2,
                 double Value_E, double Value_S,
                 double &rX1, double &rX2, double &rCOS2,
                 double &rE1, double &rE2, double &rQ,
                 double &rQElect );
};
}
///////////////////////////////////////////////////////////////////////////
#endif  // _INCLUDE_EDT1D_H_
