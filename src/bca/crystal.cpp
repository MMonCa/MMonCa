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
// CRYSTAL.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "crystal.h"

using std::endl;
using std::ostream;
using std::istream;
using std::cerr;
using std::ofstream;

namespace BCA {
char CenteredTypeStr[][15] =
     { "PRIMITIVE", "BODY CENTERED", "A END CENTERED", "B END CENTERED",
       "C END CENTERED", "FACE CENTERED" };

//-------------------------------------------------------------------------
ostream &operator<<( ostream &Out, XtalDefinition A)
{
  if(A.X==1) 
    Out << "    Atom " << A.Type << " : ";
  else
    Out << "    Atom " << A.Type << " (" << A.X*100 <<"%)/"
                   << A.Type2 << " (" << (1-A.X)*100 <<"%) : " ;
  Out << CenteredTypeStr[ (int)A.Center-1 ];
  Out << A.Origin << " "; 
  if(A.X==1)
    Out << A.BindingEnergy << " eV\n";
  else
    Out << A.BindingEnergy << "/" << A.BindingEnergy2 <<" eV\n";
  return Out;
}
//-------------------------------------------------------------------------
istream &operator>>( istream &In, XtalDefinition &A )
{
  int i;

  In >> A.Type >> A.Type2;
  In >> A.X;
  In >> i;
  A.Center = (CenteredType) i;;
  In >> A.Origin;
  In >> A.BindingEnergy >> A.BindingEnergy2;

  return In;
}
//-------------------------------------------------------------------------
int Crystal::MakeIt1( Vector Origin, CenteredType Type, int  AtomType1,
                        int  AtomType2, double X )
{
 // Variables
  int i, j, k;
  Vector V;

  // For debugging
  assert( Type>=PRIMITIVE && Type<=FACE_CENTERED );

 // Add the lattice site index
  NeighIndex++;
  if(NeighIndex > NEIGHMAX)
   {
    cerr << "# Error: Crystal::MakeIt1, Maximum number of lattice sites reached "
     << NeighIndex << " " << NEIGHMAX << endl;
    exit(MAX_LATTICE_SITES_ERROR);
   };
  int LS = NeighIndex;

 // Calcule the dimensions of the cube
  int DMax = (int) (2+2*(XTalSize+0.5));

 // Phase 1.- Build a basic cubic structure
  ncellsx = (int) DMax/2;
  ncellsy = (int) DMax/2;
  ncellsz = (int) DMax/2;

  // Put into the proper suboctant;
  V   = Origin;
  Origin  = Origin - V.Integer( TRANSLATION_PARAMETER );

  // To correct the axes orientation
  Origin.Z= -Origin.Z;
  for (k=-ncellsz; k<ncellsz; k++)
   for (j=-ncellsy; j<ncellsy; j++)
    for (i=-ncellsx; i<ncellsx; i++)
    {
      // PRIMITIVE cubic structure
       Add(Atom(AtomType1, LS, V( i,    j,        - k     ) + Origin, AtomType2,X));

     switch ( Type )
     {
      case BODY_CENTERED:
       Add(Atom(AtomType1, LS, V( i + 0.50, j + 0.50, -(k + 0.50) ) + Origin,AtomType2,X ));
       break;
      case A_END_CENTERED:
       Add(Atom(AtomType1, LS, V( i + 0.50, j + 0.50, - k     ) + Origin,AtomType2,X ));
       break;
      case B_END_CENTERED:
       Add(Atom(AtomType1, LS, V( i   , j + 0.50, -(k + 0.50) ) + Origin,AtomType2,X ));
       break;
      case C_END_CENTERED:
       Add(Atom(AtomType1, LS, V( i + 0.50, j       , -(k + 0.50) ) + Origin,AtomType2,X ));
       break;
      case FACE_CENTERED:
       Add(Atom(AtomType1, LS, V( i + 0.50, j + 0.50, - k     ) + Origin,AtomType2,X ));
       Add(Atom(AtomType1, LS, V( i   , j + 0.50, -(k + 0.50) ) + Origin,AtomType2,X ));
       Add(Atom(AtomType1, LS, V( i + 0.50, j       , -(k + 0.50) ) + Origin,AtomType2,X ));
       break;
     };
    };
  return N;
}
//-------------------------------------------------------------------------
int Crystal::MakeIt2( double a, double b, double c, double Alpha, double Beta, double Gamma )
{
  // For debugging
  assert( a>0 && b>0 && c>0 );
  assert( Alpha>0 && Beta>0 && Gamma>0 );

  // Phase 2.- Transform to the new reference system
  Vector NewX, NewY, NewZ;
  double cosA = cos( Alpha * DEGRAD );
  double cosB = cos( Beta  * DEGRAD );
  double cosC = cos( Gamma * DEGRAD );
  double sinC = sin( Gamma * DEGRAD );
  double cY = (cosA-cosB*cosC)/sinC;

   // Make the new reference vectors
  NewX( 1.0,  0.0,  0.0 );
  NewY( cosC, sinC, 0.0 );
  NewZ( cosB, cY,   sqrt( 1 - cosB*cosB - cY*cY ) );

   // Re-scale axes
  NewX = a * NewX;
  NewY = b * NewY;
  NewZ = c * NewZ;

/////
  NewX = Precision( NewX, DELTA1, 0.0 );
  NewY = Precision( NewY, DELTA1, 0.0 );
  NewZ = Precision( NewZ, DELTA1, 0.0 );

   // Calculate Cell Volume
  CellVolume = a*b*c
    *sqrt( 1 - cosA*cosA - cosB*cosB - cosC*cosC + 2*cosA*cosB*cosC );

   // Calculate Mean Atomic Radius
  int NAtoms = N / (8*ncellsx*ncellsy*ncellsz);
  MeanAtomicRadius = exp(log( 0.75*CellVolume/(M_PI*NAtoms))/3.0);
  MeanAtomicRadius *= MeanAtomicRadius; // To work with the square

   // Calculate double Density
  Dens = NAtoms / CellVolume * 1e24; // In at/cm^3

   // Conversion main loop
  LinkObject<Atom>* Item;
  Item = PointerToFirstData();
  do
   {
    Vector R = Item->Data.R;
    R = ChangeToReference( NewX, NewY, NewZ, R );
    R = Precision( R, DELTA1, 0.0 );
    Item->Data.R = R;
    Item = PointerToNextData();
   }
  while( Item != NULL );
  NeighIndex = -1; // Second initialization for Neigh Index
  return( N );
}
//-------------------------------------------------------------------------
void Crystal::CreateEDT( char *EDTFile, EDT1D Dens1D[N_ATOMS+1] )
{
  double I,J,K;
  LinkObject<Atom> * Item;
  time_t T0,T1;
  
  time(&T0);
  
  double LX = LatticeParameter.X;  // In (A)
  double LY = LatticeParameter.Y;  //
  double LZ = LatticeParameter.Z;  //
  
  DENS.InitEmptyTable( LX,LY,LZ,NX,NY,NZ);

  int NN=NX;
  for( I=-LX/2; I<LX/2; I+=LX/NX )
  {
   for( J=-LY/2; J<LY/2; J+=LY/NY )
    for( K=-LZ/2; K<LZ/2; K+=LZ/NZ )
     {
      double Value = 0; 
      
      Item=PointerToFirstData();
      do
        {
         double X = (I-Item->Data.R.X); // In (A)
         double Y = (J-Item->Data.R.Y); 
         double Z = (K-Item->Data.R.Z);     
         double R = sqrt( X*X + Y*Y + Z*Z );   

         // Weighted Mean value of the densities when various atoms are
         //  at one lattice site
         if(Item->Data.X==1)
          Value+=Dens1D[ Item->Data.Index ].Interpolate1D( R );
         else
          Value+=Dens1D[ Item->Data.Index1 ].Interpolate1D( R )
                    *Item->Data.X +
                 Dens1D[ Item->Data.Index2 ].Interpolate1D( R )
                    *(1-Item->Data.X);
         Item = PointerToNextData();
        }
      while( Item != NULL );
        
      DENS.SetValue( (int)((I+LX/2)*NX/LX), 
                     (int)((J+LY/2)*NY/LY), 
                     (int)((K+LZ/2)*NZ/LZ), Value );    
     };
    cerr.width(8);
    cerr << NN-- << "\b\b\b\b\b\b\b\b\b"; 
  };

  char *EDTDIR, NewFileName2[255];
  EDTDIR = getenv( "EDTDIR" ); // Read the environment variable
  if( EDTDIR != NULL ) strcpy( NewFileName2, EDTDIR );
  else                 strcpy( NewFileName2, EDTDIRDEFAULT );
  strcat( NewFileName2, EDTFile );

  time(&T1);
  cerr << "# "<< difftime(T1,T0) << " seconds\n";
  
//  DENS.SaveASCIITable3D( NewFileName2 );
  DENS.SaveBinTable3D( NewFileName2 );
}
//-------------------------------------------------------------------------
int Crystal::Update( double XS, double IR, Vector LP,
             double GL, double SD, double IP)
{
  LatticeParameter     = LP;        // Measured in Angstroms
  if( LatticeParameter == NullVector ) LatticeParameter(5.431,5.431,5.431);
  XTalSize             = XS;        //      in lattice units
  InteractionRadius    = IR;        //      in A
  GhiLimit             = GL;        //      in lattice units
  SimultaneousDistance = SD;        //      in A
  IR2                  = InteractionRadius*InteractionRadius;
  MaxIP                = IR2;

  DeleteAll();
  Set( (List<Atom>) Neighbors[0] );

  return(1);
}
//-------------------------------------------------------------------------
// Make Neighbors /////////////////////////////////////////////////////////
int Crystal::MakeNeighbors( Vector Origin /* Lattice units */ )
{
  NeighIndex++;
  if(NeighIndex > NEIGHMAX)
   {
    cerr << "# Error: Crystal::MakeNeighbors, Maximum number of lattice sites reached\n";
    exit(MAX_LATTICE_SITES_ERROR);
   };
  Origin( Origin.X*LatticeParameter.X,
      Origin.Y*LatticeParameter.Y,
     -Origin.Z*LatticeParameter.Z);  // Correct the axis orient.

  double IR  = InteractionRadius+( 0.5*sqrt(3.0)) * LatticeParameter.GetMax();
  LinkObject<Atom>* Auxiliar = PointerToFirstData();
  do
  {
   if( (double )(Auxiliar->Data.R - Origin) < (IR) )
    {
      Auxiliar->Data.R = Auxiliar->Data.R - Origin;
      Neighbors[NeighIndex].Add( Auxiliar->Data );
    };
   Auxiliar = PointerToNextData();
  } while (Auxiliar != NULL);

  int i = Neighbors[NeighIndex].Sort( ByDistance, NullVector, sortASCENDENT );
  return i;
}
//-------------------------------------------------------------------------
int Crystal::Rotate( double Angle, AxeType Axe, Vector Reference = NullVector)
{
  // Variables
  Vector AxeX, AxeY, AxeZ;

  if(Angle==0.0) return(N);

  // For debugging
  assert( (Axe==AXE_X)||(Axe==AXE_Y)||(Axe==AXE_Z) );

  // Selection of AXEs for Rotation
  switch (Axe)  // Care!! This triedro is not direct
   {
    case AXE_X: AxeX( 1.0,         0.0,        0.0 );
                AxeY( 0.0,  cos(Angle), sin(Angle) );
                AxeZ( 0.0, -sin(Angle), cos(Angle) );
                break;
    case AXE_Y: AxeX( cos(Angle), 0.0, -sin(Angle) );
                AxeY(        0.0, 1.0,         0.0 );
                AxeZ( sin(Angle), 0.0,  cos(Angle) );
                break;
    case AXE_Z: AxeX(  cos(Angle), sin(Angle), 0.0 );
                AxeX( -sin(Angle), cos(Angle), 0.0 );
                AxeZ(         0.0,        0.0, 1.0 );
                break;
   };

  // It Does the rotation. OPTIMIZED
  int Index = 1;
  LinkObject<Atom>* BestData = PointerToFirstData();
  do
   {
    BestData->Data.R =  ChangeToReference(
                       AxeX, AxeY, AxeZ,
                       BestData->Data.R - Reference
                     )
            + Reference;

    Index++;
    BestData = PointerToNextData();
   }
  while ( Index <= N );
  return( Index );
}
//-------------------------------------------------------------------------
// Rotate crystal (multiply by rotation matrix) OPTIMIZED /////////////////
//  Select the Neighbor list gived by Near.LatticeSite
int Crystal::Rotate( double M[3][3] )
{
  int Index;
  Vector R ;
  LinkObject<Atom>* Item,* Item2;
  List<Atom> *Old = &Neighbors[ 0  ];
 // List<Atom> *Old = &Neighbors[ Near.LS ];

  Item    = Old->PointerToFirstData();
  Item2   = PointerToFirstData();               // First element
  for( Index=0; Index < Old->N; Index++ )   
   {
    R = Item->Data.R;                   // Coordinates
    R(                              // Rotate coordinates
      R.X*M[0][0] + R.Y*M[0][1] + R.Z*M[0][2],
      R.X*M[1][0] + R.Y*M[1][1] + R.Z*M[1][2],
      R.X*M[2][0] + R.Y*M[2][1] + R.Z*M[2][2]
     );

    Item2->Data.R =  R;                 // Update
    Item    = Old->PointerToNextData();
    Item2   = PointerToNextData();
   };
 
  return( Index );
}
//-------------------------------------------------------------------------
int Crystal::UnRotate( )
{
  int Index;
  Vector R ;
  LinkObject<Atom>* Item,* Item2;
  List<Atom> *Old = &Neighbors[ 0  ];
  Item    = Old->PointerToFirstData();
  Item2   = PointerToFirstData();               // First element
  for( Index=0; Index < Old->N; Index++ )   
   {
    R = Item->Data.R;                   // Coordinates
    Item2->Data.R =  R;                 // Update
    Item    = Old->PointerToNextData();
    Item2   = PointerToNextData();
   }; 
  return( Index );
}
//-------------------------------------------------------------------------
// Print the crystallite in XYZ format to represent with XMOL /////////////
int Crystal::ShowToXMol(const std::string &FileName, Vector VP )
{
  int i = 0;
  ofstream Out;
  LinkObject<Atom>* Item;

  Out.open( FileName.c_str(), std::ios::out );
  if( !Out ) 
   { cerr << "# Error: Crystal::ShowToXMol, Can not open the file " 
          << FileName << endl; 
     exit(READ_ERROR);
   }
  else
  {
   Out << N << endl << endl;
   Item = PointerToFirstData();
   do
   {
    AtomDefinition A = Atoms.Search( Item->Data.Index )->Data;
    Out << A.GetName();
    Item->Data.PrintXYZ( Out );
    Item = PointerToNextData();
    i++;
   } while(Item != NULL);
//   Out << VP.X << " " << VP.Y << " " << VP.Z << endl;
   Out.close();
  };
  return( i );
}
//-------------------------------------------------------------------------
ostream &operator<<( ostream &Out, Crystal X )
{
  Out << (List<Atom>) X;
  return( Out );
}
//-------------------------------------------------------------------------
// Print a List of ATOMS in XYZ format to represent with XMOL /////////////
int ShowToXMol( List<Atom> LIST ,char* FileName )
{
 LinkObject<Atom>* Item;
 ofstream Out;

 Out.open( FileName, std::ios::out );
 if( !Out ) 
 { cerr << "# Error: ShowToXMol, Can't open the file " 
        << FileName << endl; 
   exit(READ_ERROR);
 }
 else
 {
  Out <<  LIST.N << endl << endl;
  Item = LIST.PointerToFirstData();
  do
  {
    AtomDefinition A = Atoms.Search( Item->Data.Index )->Data;
    Out << A.GetName()<< "3";
    Item->Data.PrintXYZ( Out );
    Item = LIST.PointerToNextData();
  } while(Item != NULL);
  Out << endl;
 };
 Out.close();
 return( LIST.N );
}
//-------------------------------------------------------------------------
// Print a List of PROJECTILES in XYZ format to represent with XMOL ///////
int ShowToXMol( List<Projectile> LIST ,char* FileName, int Option )
{
 LinkObject<Projectile>* Item;
 ofstream Out;

 Out.open( FileName, std::ios::out );
 if( !Out ) 
 { cerr << "# Error: ShowToXMol; Can't open the file " 
        << FileName << endl; 
   exit(READ_ERROR);
 }
 else
 {
  Out <<  (LIST.N * 2 - 1) << endl << endl;
  Item = LIST.PointerToFirstData();
  do
  {
    Item->Data.PrintXYZ( Out, Option );
    Item = LIST.PointerToNextData();
  } while(Item != NULL);
  Out << endl;
 };
 Out.close();
 return( LIST.N );
}
//-------------------------------------------------------------------------
int Crystal::SearchTargets( List<Projectile> &Targets,
                Projectile &CurrentProjectile,
                double       &IP2Min, // Lower Impact Parameter
                Random     &RND,
                double       FrontSurface,
                double       BackSurface,
                List<Atom> &Lista2 
              )
{
  // Variables
  Vector DeltaR, Org, CurrentDir, CurrentR;

  LinkObject<Projectile>* PossibleTarget2;
  LinkObject<Atom>* PossibleTarget;

  double   MinGhi = 4096.0;
  double   MaxGhi = 4096.0;
  double   Ghi;
  double   OldGhi = CurrentProjectile.Ghi + DELTA1;
  double   IP2Near    = 4096.0;
         IP2Min     = 4096.0; // apac(7)
  double   IP2;

  double   OutsideGhi = 4096.0;

  Atom   NearOut;

  // Code
  assert( OldGhi>=GhiLimit );

  CurrentDir = CurrentProjectile.Dir;
  CurrentR   = CurrentProjectile.R;
  Org        = CurrentProjectile.Abs - CurrentR;

  // I delete nearest atom list
    LP.DeleteAll();

  // First pass -----------------------------------------------------------

  PossibleTarget = PointerToFirstData();
  do
   {
      // Get a copy of atom to modify
      Projectile Item = PossibleTarget->Data;

      Item.InitPos = Item.R;
      
      // Select the atom in these lattice site (Al_{x}Ga_{1-x}As -like composites
      if( RND.Generate() < Item.X ) Item.Index = Item.Index1;
      else                          Item.Index = Item.Index2;

      // Thermal displacemet
      Item.ThermalDisplacement( Item.Index, RND);

    // TEST WHETHER THE PARTNER IS IN THE CORRECT FORWARD RANGE
    DeltaR = Item.R - CurrentR;

    // Add absolute positions within the interaction sphere
    if( (double )DeltaR <= RI ) LP.Add( Org+Item.R ); 

    Ghi    = CurrentDir * DeltaR;
    if( ( Ghi > OldGhi ) && (Ghi <= MaxGhi ) )
     {
      // CHECK IMPACT PARAMETER
      IP2 = DeltaR*DeltaR - Ghi*Ghi;
      if( IP2 < MaxIP )
       {
    // Fourth condition (X >= 0) Begin of CRYSTAL
    double XPosition = OR.ProyectaX( Org + Item.R );
    if (
         ( XPosition >= FrontSurface ) && ( XPosition <= BackSurface )
       )
     {

      // FIND THE FIRST Significant encounter along the projectile track
      if( Ghi < MinGhi )
       {
        MinGhi  = Ghi;
        MaxGhi  = Ghi + SimultaneousDistance;
        IP2Near = IP2;
        Near    = Item;
       };
      Targets.Add( Item );
     };
       }
      // NO COLLISION: FIND MINIMUM IMPACT PARAMETER FOR OUTSIDE PARTNERS
      else
       {
     if( IP2 < IP2Min ) // For outside partners
      {
        IP2Min     = IP2;   // Store minimum impact parameter
        OutsideGhi = Ghi;   //  and Ghi parameter
        NearOut    = Item;
      };
       };
     };
    PossibleTarget = PointerToNextData();
   }
  while ( PossibleTarget != NULL );

// This section is deactivated not to take into account 
//  interstitials in a deterministically manner
#ifdef INTERSTICIALES
// Search in interstitials list -------
if(Lista2.N!=0)
{
  PossibleTarget = Lista2.PointerToFirstData();
  do
   {
      // Get a copy of atom to modify
      Projectile Item = PossibleTarget->Data;

      Item.InitPos = Item.R;

      // Thermal displacemet
      Item.ThermalDisplacement( Item.Index, RND);

    // TEST WHETHER THE PARTNER IS IN THE CORRECT FORWARD RANGE
    DeltaR = Item.R - CurrentR;

    // Add absolute positions within the interaction sphere
    if( (double )DeltaR <= RI ) LP.Add( Org+Item.R ); 

    Ghi    = CurrentDir * DeltaR;
    if( ( Ghi > OldGhi ) && (Ghi <= MaxGhi ) )
     {
      // CHECK IMPACT PARAMETER
      IP2 = DeltaR*DeltaR - Ghi*Ghi;
      if( IP2 < MaxIP )
       {
    // Fourth condition (X >= 0) Begin of CRYSTAL
    double XPosition = OR.ProyectaX( Org + Item.R );
    if (
         ( XPosition >= FrontSurface ) && ( XPosition <= BackSurface )
       )
     {

      // FIND THE FIRST Significant encounter along the projectile track
      if( Ghi < MinGhi )
       {
        MinGhi  = Ghi;
        MaxGhi  = Ghi + SimultaneousDistance;
        IP2Near = IP2;
        Near    = Item;
       };
      Targets.Add( Item );
     };
       }
      // NO COLLISION: FIND MINIMUM IMPACT PARAMETER FOR OUTSIDE PARTNERS
      else
       {
     if( IP2 < IP2Min ) // For outside partners
      {
        IP2Min     = IP2;   // Store minimum impact parameter
        OutsideGhi = Ghi;   //  and Ghi parameter
        NearOut    = Item;
      };
       };
     };
    PossibleTarget = Lista2.PointerToNextData();
   }
  while ( PossibleTarget != NULL );
}
#endif

  // Second pass -----------------------------------------------------------

  if( Targets.N > 0 )            // There are targets //
  {

    IP2Min = IP2Near;

    PossibleTarget2 = Targets.PointerToFirstData();
    do
    {
      DeltaR = PossibleTarget2->Data.R - CurrentR;
      Ghi    = CurrentDir * DeltaR;
      IP2    = DeltaR*DeltaR - Ghi*Ghi;

      if( Ghi > MinGhi )
      {
     if( Ghi <= MaxGhi )
     {
       // Delete targets outside interaction sphere

        double DQ = IR2 - ( Ghi-MinGhi ) * ( Ghi-MinGhi );
        if ( (IP2>DQ) || (IP2Near>DQ) )
        {
          Targets.Delete(PossibleTarget2);
          continue;
        };
     }
     else
     {
       // Delete targets further MaxGhi

       Targets.Delete(PossibleTarget2);
       continue;
     };
      };

      // Solve the Minimum Impact Parameter
      if( IP2 < IP2Min ) IP2Min = IP2;

    }
    while (
       PossibleTarget2 = Targets.PointerToNextData(),
       PossibleTarget2 != NULL
      );
  }
  else                       // There are no targets //
  {
    if (OutsideGhi == 4096.0 ) OutsideGhi = GhiLimit;
    CurrentProjectile.Ghi = OutsideGhi; // Update Ghi
    Near = NearOut;
  };

  return ( Targets.N );
}
}
//------------------------------------------------------------------------- 
