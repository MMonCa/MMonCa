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
// EDT1D.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "edt1d.h"

using namespace std;
namespace BCA {
//-------------------------------------------------------------------------
EDT1D::EDT1D()
{
 Radius[0] = 0;
 Density[0] = 0;
 ND = 1;
}
//-------------------------------------------------------------------------
void EDT1D::ReadTable(const std::string &FileName, int PP )
{
    int i;
    ifstream FP;
    char *EDTDIR, NewFileName[255];

    if( PP ) EDTDIR = getenv( "EDTDIR" );
    if( EDTDIR != NULL ) strcpy( NewFileName, EDTDIR );
    else                 strcpy( NewFileName, EDTDIRDEFAULT );
    strcat( NewFileName, FileName.c_str());

#ifdef WINDOWS
    FP.open(NewFileName,ios::in|ios::nocreate);
#else
    FP.open(NewFileName,ios::in);
#endif
    if (!FP){
      cerr << endl
           << "# Error: EDT1D::ReadTable, File (1D) not found: "
           << NewFileName << endl;
      exit( READ_ERROR );
     };

    i=0;
    while (FP.eof()==0)
    {
        FP >> Radius[i] >> Density[i];       
        if(Radius[i]<0.0 || Density[i]<=0) break;
        Density[i] = log(Density[i]);  // In logarithmic scale
        i++;
        if (i==NMXV) break;
    }
    ND=i-1;
    if( PP ) cerr << "(" << NewFileName << ") ND = " << ND << endl;

    FP.close();
}
//-------------------------------------------------------------------------
double EDT1D::GetMaxR() { return Radius[ND]; }
//-------------------------------------------------------------------------
double EDT1D::Interpolate1D( double r )  // r measured in (A)
{
    int i;
    double x0,x1,y0,y1;

    if(r<=0) return 0;

    // Sequencial search !! //for(i=1;i<=ND;i++) if (Radius[i]>r) break;
    // Binary search !!
    int lo=1;
    int hi=ND;
    do{
        i = ( hi + lo )/ 2;
        if (Radius[i]>r) hi = i;
        else             lo = i;
      }while (hi-lo != 1);
    i=hi;

  #ifndef NOEXTRAPOLATION    // Solve WITH    extrapolation
    if (i>=ND) i=ND;
  #else                      // Solve WITHOUT extrapolation
    if (i>=ND) return (1e-10);
  #endif
    x0=Radius [i-1];
    x1=Radius [i  ];
    y0=Density[i-1];
    y1=Density[i  ];
    y0=exp((r-x0)*(y1-y0)/(x1-x0)+y0);
    return( y0 ); // in atomic units
}
//-------------------------------------------------------------------------
double EDT1D::GetRs( double R )
{
  double Value = Interpolate1D( R );
  assert( Value != 0.0 );
  return( pow( 0.75/M_PI/Value,1./3.)/bohr_radius*1e-10 );
}
}
//-------------------------------------------------------------------------
