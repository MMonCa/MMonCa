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
// LAYER.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "layer.h"

using std::ostream;
using std::endl;
using std::cout;
using std::cerr;

namespace BCA {
// Minimum Common Multiply
int MCM( int a, int b )
{
 int m = 1;
 if( b==0 ) return a;
 for( int i=a>b?a:b;i>1;i--)
  if( a%i==0&&b%i==0) { m=m*i; a=a/i;b=b/i;}
 return m;
} 

// Output stream overloaded
ostream &operator<<( ostream &Out, LayerDefinition A)
{
  Out << "\n==> ";
  switch (A.Amorphous)
   {
    case 0 : Out << "CRISTALLINE layer";      break;
    case 1 : Out << "POLYCRISTALLINE layer";  break;
    case 2 : Out << "AMORPHOUS layer";        break;
    default: Out << "UNKNOWN layer";
   };
  Out << " ( " << A.XMin << " to " << A.XMax << " ) A" << endl;
  Out << "    LatticeParameter " << A.LatticeParameter << " A" << endl;
  Out << "    Angles           " << A.Angles    << endl;
  Out << "    EDTFile          ' " << A.EDTFile << "'";
  if( A.EDTCreate ) Out << " (new)\n"; else Out <<" (old)\n";
  Out << A.XTal;
  return Out;
}

///////////////////////////////////////////////////////////////////////////
// Class Layer
// Make all layers in simulation ----------------------------------

void LayerClass::MakeLayers( List<LayerDefinition> &Bulk, int Initialize, 
                                int Show )
{
 LinkObject<LayerDefinition>* Item;
 LinkObject<XtalDefinition>*  Item2;

 List<XtalDefinition>         XTal;
 Crystal *            Crystallite;
 Vector               LatticeParameter;
 Vector               Angles;
 int i, nLayer = 0;

 // Default
 if( Bulk.NumberOfElements() == 0 ) // SILICON by Default
 {
  Vector     Origin;
  LayerDefinition L; // Default values are for Silicon
  XtalDefinition  X;

  Origin(0.00,0.00,0.00);
  X.Set( Origin, 1, FACE_CENTERED );
  L.XTal.Add( X );

  Origin(0.25,0.25,0.25);
  X.Set( Origin, 1, FACE_CENTERED );
  L.XTal.Add( X );

  Bulk.Add( L );
 };

 FrontSurface =  1.0e10;
 BackSurface  = -1.0e10;

 Item = Bulk.PointerToFirstData();
 // Main loop
 do
  {
   Crystallite = new Crystal( XTalSize, InteractionRadius, NullVector,
                   GhiLimit, SimultaneousDistance );

   // Get the minimum and maximum surfaces
   if( Item->Data.XMin < FrontSurface ) FrontSurface = Item->Data.XMin;
   if( Item->Data.XMax > BackSurface  ) BackSurface  = Item->Data.XMax;

   // Some variables of each layer
   XTal             = Item->Data.XTal;
   LatticeParameter = Item->Data.LatticeParameter;
   Angles           = Item->Data.Angles;
   // Do the xtal
   Item2 = XTal.PointerToFirstData();
   do
   {
       Crystallite->MakeIt1( Item2->Data );
       Item2 = XTal.PointerToNextData();
   }
   while( Item2 != NULL );
   // Xtal transformation
   Crystallite->MakeIt2(
              LatticeParameter.X,
              LatticeParameter.Y,
              LatticeParameter.Z,
              Angles.X,
              Angles.Y,
              Angles.Z
              );
   // Outputs xtal information (if requested)
   nLayer++;
   char sss[200];
   sprintf(sss,"Layer_%02d.xyz",nLayer);
   if(Show)
   {
    cout << Item->Data;
    cout << "    Cell volume        = " << Crystallite->CellVolume << " A^3" << endl;
    cout << "    Mean atomic radius = " << sqrt(Crystallite->MeanAtomicRadius) << " A"<< endl;
    cout << "    Theorical density  = " << Crystallite->Dens << " at/cm^3" << endl;
    cout << endl;
    Crystallite->ShowToXMol( sss );   // Show crystallite in .XYZ format 
   }
  
   // Generate the list of neighbors
   Item2 = XTal.PointerToFirstData();
   do
   {
     Crystallite->MakeNeighbors( Item2->Data );
     Item2 = XTal.PointerToNextData();
   }
   while( Item2 != NULL );

  // Updates the crystallite and set some variables
  Crystallite->Update( XTalSize, InteractionRadius, LatticeParameter,
              GhiLimit, SimultaneousDistance );

  // Create the IADS electron density
  static EDT1D Dens1D[N_ATOMS+1];
  int NFN = (Item->Data.EDTFile[0]=='\0');
    
  char N[N_ATOMS+1], PP[2];
  if( NFN ) strcpy( Item->Data.EDTFile, "EDT_" );

  for( i=0;i<N_ATOMS; i++) N[i] = 0;
  Item2 = XTal.PointerToFirstData();
  do 
   { 
    N[ Item2->Data.Type ] ++; 
    N[ Item2->Data.Type2 ] ++; 
    Item2 = XTal.PointerToNextData(); 
   }   
  while( Item2 != NULL );
  int mcm = N[1]; for( i=2; i<N_ATOMS; i++ ) mcm = MCM( mcm, N[i] );

  for( i=1; i<N_ATOMS; i++ )
   if( N[i]>0 ) 
    {
     if( Item->Data.EDTCreate!=0)
      {
       char FN[30];
       sprintf( FN, "N%02i.den", Atoms.Search(i)->Data.Z ); 
       cerr << "# Reading (" << FN << ") ";
       Dens1D[i].ReadTable( FN );
      }   
     if( NFN )   
      {
       strcat( Item->Data.EDTFile, Atoms.Search(i)->Data.GetName());
       PP[0]=N[i]/mcm+48;          PP[1]=0;
       if( PP[0] > '1' ) strcat( Item->Data.EDTFile, PP );
      }; 
    };
  
  if( !NFN ) Item->Data.EDTFile[ strlen(Item->Data.EDTFile) - 1 ] = '\0';
  if( Item->Data.EDTCreate!=0 )
   {
    cerr << "# Creating IADS in " << Item->Data.EDTFile << endl;
    Crystallite->CreateEDT( Item->Data.EDTFile, Dens1D );
   }
  else 
   {
    if( Initialize ) 
     if (Crystallite->DENS.ReadTable3D( Item->Data.EDTFile, 1 ))
      {
       // Create 3D table density
       for( i=1; i<N_ATOMS; i++ )
        if( N[i]>0 ) 
         {
           char FN[30];
           sprintf( FN, "N%02i.den", Atoms.Search(i)->Data.Z ); 
           cerr << "# Reading (" << FN << ") ";
           Dens1D[i].ReadTable( FN );       
         };
       cerr << "# Creating Isolate Atom Density Superposition in " << Item->Data.EDTFile << endl;
       Crystallite->CreateEDT( Item->Data.EDTFile, Dens1D ); 
      };  

   }; 

  // ADD this layer
  Layers.Add( *Crystallite );

  // Next ...
  Item = Bulk.PointerToNextData();
 } while ( Item!=NULL );
 if(Show)
 {
  cout << "    Front surface at " <<  FrontSurface << " A" << endl;
  cout << "    Back surface at " << BackSurface  << " A" << endl;
  cout << endl;
 }
}
}
///////////////////////////////////////////////////////////////////////////
