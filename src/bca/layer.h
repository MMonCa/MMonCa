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
// LAYER.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_LAYER_H_
#define _INCLUDE_LAYER_H_

#include "vector.h"
#include "list.h"
#include "crystal.h"
namespace BCA {
int MCM( int a, int b );

// Class Layer definition /////////////////////////////////////////////////

class LayerDefinition
{
 public:
   Vector   LatticeParameter;
   Vector   Angles;
   int      Amorphous;
   double     XMin;
   double     XMax;
   char         EDTFile[255];
   int          EDTCreate;  // (=0) Don't create, read it. (!=0) Create and save it
   List<XtalDefinition> XTal;

   LayerDefinition()    ////////// CONSTRUCTOR
   {
     LatticeParameter(  5.431, 5.431, 5.431 );
     Angles      ( 90.0,  90.0,  90.0   );
     Amorphous   =   0;
     XMin    =   0.0;   // (A)
     XMax    = 1e9;     // (A)
     EDTFile[0]  = '\0';
     EDTCreate   = 0; // Don't create file
   };

   void Default()   // Reset the current object
   {
     List<XtalDefinition> NewXTal;
     LatticeParameter(  5.431, 5.431, 5.431 );
     Angles      ( 90.0,  90.0,  90.0   );
     Amorphous   =   0;
     XMin    = 0.0;
     XMax    = 1e9;
     XTal    = NewXTal;
     EDTFile[0]  = '\0';
     EDTCreate   = 0; // Don't create file
   };

   friend std::ostream &operator<<(std::ostream &Out, LayerDefinition A);
};

///////////////////////////////////////////////////////////////////////////
//
// Class Layer
//
class LayerClass
{
    double       XTalSize;
    double       InteractionRadius;
    double       GhiLimit;
    double       SimultaneousDistance;
  public:
    double       BackSurface;
    double       FrontSurface;
    List<Crystal> Layers;

 void Init( double XS=0.5, double IR=2.7155, double GL=0.1, double SD=0.5 )
 {
   XTalSize    = XS;        //  in lattice units
   InteractionRadius = IR;      //  in (A)
   GhiLimit    = GL;        //  in lattice units
   SimultaneousDistance = SD;       //  in (A)
 };
 // Make all layers in simulation ----------------------------------
 void MakeLayers( List<LayerDefinition> &Bulk, int Initialize = 1, 
                    int Show =0 );

}; // End of LayerCLASS
}
///////////////////////////////////////////////////////////////////////////
#endif  // _INCLUDE_LAYER_H_
