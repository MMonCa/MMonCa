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
// VECTOR.CPP
//
///////////////////////////////////////////////////////////////////////////

#include "vector.h"

using std::endl;
using std::cout;
using std::ostream;
using std::istream;

namespace BCA {

// Output stream overloaded
ostream &operator<<( ostream &Out, Vector V )
{
  Out << "( ";
  Out.width(WIDTH);
  Out << V.X << ", ";
  Out.width(WIDTH);
  Out << V.Y << ", ";
  Out.width(WIDTH);
  Out << V.Z << " )";

  return Out;
}

// Input stream overloaded
istream &operator>>( istream &In, Vector &V )
{
  In >> V.X;
  In >> V.Y;
  In >> V.Z;
  return In;
}

///////////////////////////////////////////////////////////////////////////
//
//  Orientacion class
//     I define within this class a virtual coordinate system in order to
//   define crystalographic orientations different from <100>
//     If I want <110> orientation, the working coordinate system does not
//   change, avoiding 3D problems and speeding-up the code
//     I need to take into account:
//           
//       -The surface problem, it will be defined by planes
//       -Obtain the final coordinates

void Orientacion::Orienta( Vector V, Vector Flat, double Tha_Cut, double Phi_Cut, 
                            int Show )
{
// V and Flat must be perpendicular
    EjeX = Unitary( V ); 
    EjeY = Unitary( Flat );
    EjeZ = EjeX ^ EjeY;
    
    Phi_Cut = Phi_Cut * M_PI/180;
    Tha_Cut = Tha_Cut * M_PI/180;

    double cf = cos(Phi_Cut);
    double sf = sin(Phi_Cut);
    double ct = cos(Tha_Cut);
    double st = sin(Tha_Cut);
    Vector nEjeX,nEjeY,nEjeZ;
    nEjeX =   cf*ct * EjeX + sf*ct * EjeY + st * EjeZ;
    nEjeY = - sf    * EjeX + cf    * EjeY;
    nEjeZ = - st*cf * EjeX - st*sf * EjeY + ct * EjeZ;

    EjeX =   nEjeX;
    EjeY =   nEjeY; 
    EjeZ = - nEjeZ; // Indirect triedrum
    if(Show)
    {
       cout << endl;
       cout << " Axe X " << EjeX << endl;
       cout << " Axe Y " << EjeY << endl;
       cout << " Axe Z " << EjeZ << endl;
    }
}

double Orientacion::ProyectaX( Vector V )
{
    return EjeX*V;
}

double Orientacion::ProyectaY( Vector V )
{
    return EjeY*V;
}

double Orientacion::ProyectaZ( Vector V )
{
    return EjeZ*V;
}

Vector Orientacion::Proyecta( Vector V )
{
    Vector T;
    T( EjeX*V, EjeY*V, EjeZ*V );
    return T;
}

// Global variables -------------------------------------------------
Vector NullVector;  // Null vector
Vector XVector;     // X direction vector
Orientacion OR;     // Controls the orientation

// Global variables inicialization
void Initialize()
{
  NullVector(0,0,0);
  XVector(1,0,0);
}
}
///////////////////////////////////////////////////////////////////////////
