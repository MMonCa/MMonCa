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
// EDT3D.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_EDT3D_H_
#define _INCLUDE_EDT3D_H_

// #define NOEXTRAPOLATION
//  If not defined makes a extrapolation else no extrapolation
//   in Interpolate3D() and GetRs()

#include <string.h>
#include <math.h>

#include "defs.h"
#include "ssf.h"
#include "stopping.h"
#include "list.h"

namespace BCA {
extern class STP STPO;
extern List<Vector> LP;

// Class: Electronic Density Table
class EDT
{
  protected:    // EDT 3D
  
    double LADOX, LADOY, LADOZ;      // (=5.431*3 A)
    int    N_POINTSX, N_POINTSY, N_POINTSZ;  // (=112*3)
    double *DATA;       // Dynamic data
    long   SizeOfDATA;
    int    dX,dY;

  public:   // EDT 3D
    int    ReadTable3D(const std::string &FileName, int SaveIt );
    void   InitEmptyTable( double LX, double LY, double LZ,
                int nX, int nY, int nZ );
    void   SetValue( int i, int j, int k, double Value );
    void   SaveASCIITable3D( char * FileName );
    void   SaveBinTable3D( char * FileName );
    double Interpolate3D( double X, double Y, double Z );
    double NoInterpolate3D( double X, double Y, double Z );
    double GetRs3D( double X, double Y, double Z );
    void   Cut3D( double x0, double y0, double z0,
              double A,  double B,  double C  );
};

#define NPASOS_S 30

// Density CLASS

class DensityClass : public EDT
{
 public:

 double LossesNew( double Z1, double M1, Vector PAbs, Vector Dir,
           double ENERGY, double E1, double AllowedDistance, int Z2 );
};
}
#endif //_INCLUDE_EDT3D_H_
