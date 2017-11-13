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
//  TABLE3D.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "table3d.h"

using std::cerr;
using std::endl;
using std::ofstream;
using std::ifstream;

namespace BCA {

static TC MDSolver;

TBL::TBL()
{
 int i1, i2;

 for( i1=0; i1<N_ATOMS;  i1++ )
  for( i2=0; i2<N_ATOMS;  i2++ )
   {
     GRAN[ i1 ][ i2 ] = NULL;
   };
 for( i1=0; i1<HEADERSIZE; i1++ ) Cabecera[i1]='%';
 Cabecera[i1]='\0';
}

// Read Tables From Disk ///////////////////////////////////////////////////

void TBL::ReadTable( )
{
 int i_Atom, j_Atom, Z1, Z2;
 LinkObject<AtomDefinition> *Atom1,*Atom2;
 char FileName[] = "eel____.dat";

 assert( &Atoms != NULL ); // For debugging

// cout << endl;
 for( i_Atom = 0; i_Atom < Atoms.NumberOfElements(); i_Atom++ )
 {
   Atom1 = Atoms.Search( i_Atom+1 );
   Z1    = Atom1->Data.Z;
   for( j_Atom = 0; j_Atom < Atoms.NumberOfElements(); j_Atom++ )
   {
    Atom2 = Atoms.Search( j_Atom+1 );
    Z2     = Atom2->Data.Z;
    GetSSFName(Z1,Z2,FileName);
    if( !LoadBINARY(i_Atom,j_Atom,FileName) )
       if( LoadASCII(i_Atom,j_Atom,FileName)) 
           SaveBINARY(i_Atom,j_Atom,FileName);
   }    
 }
// cout << endl;
}

///////////////////////////////////////////////////////////////////////////

int TBL::LoadBINARY( int nZ1, int nZ2, char * FileName )
{
 ifstream INPUT;
 char NewFileName[200];
 union { short int w; char b; } endian;
 endian.w = 0x0100;
 
 if(endian.b==1)    // BIG Endian
  {
      FileName[ 8] = 'b';
      FileName[ 9] = 'i';
      FileName[10] = 'g';   // Trying HP/UX or SUN   architecture
  }
 else
  {
      FileName[ 8] = 'l';
      FileName[ 9] = 'i';
      FileName[10] = 't';   // Trying ALPHA or INTEL architecture
  };
  TABLEDIR = getenv( "TABLES3DDIR" );  
  if( TABLEDIR != NULL ) strcpy( NewFileName, TABLEDIR );
  else                   strcpy( NewFileName, TABLEDIRDEFAULT );
  strcat( NewFileName, FileName );
#ifdef WINDOWS
  INPUT.open( NewFileName, ios::binary|ios::nocreate );
#else
  INPUT.open( NewFileName ); 
#endif
  if( INPUT )
     {
       tTABLE * pTABLE;

       pTABLE = new tTABLE ;
       if( !pTABLE )
       {
        cerr << endl;
        cerr << "# Error: TBL::LoadBinary, I Can't Allocate " << sizeof(tTABLE)
         << " bytes !!" << endl;
        exit(MEMORY_ERROR);
       };

#ifdef GCC3
       INPUT.read( (char *) Cabecera, HEADERSIZE );
       INPUT.read( (char *) pTABLE->TABLE, sizeof(tTABLE) );
#else      
       INPUT.read( (unsigned char *) Cabecera, HEADERSIZE );
       INPUT.read( (unsigned char *) pTABLE->TABLE, sizeof(tTABLE) );
#endif
       GRAN[nZ1][nZ2] = pTABLE;
       INPUT.close();
       
       cerr << "# Found " << _BLUE_ << NewFileName << " " << _NORMAL_
            << "(" << INPUT.gcount()/sizeof(tTABLE)*100 << "%) OK. ";

       //cerr << Cabecera << endl;
       cerr << endl;       

       return 1;
    }
    else return 0;
}

///////////////////////////////////////////////////////////////////////////

int TBL::LoadASCII( int nZ1, int nZ2, char *FileName )
{
  char NewFileName[255];
  ifstream INPUT;
  int index_E, index_S, index_PHI, index_THETA, index_ALPHA, index_V;      
  double Value;
 
  FileName[ 8] = 'd';
  FileName[ 9] = 'a';
  FileName[10] = 't';

  TABLEDIR = getenv( "TABLES3DDIR" );  
  if( TABLEDIR != NULL ) strcpy( NewFileName, TABLEDIR );
  else                   strcpy( NewFileName, TABLEDIRDEFAULT );
  strcat( NewFileName, FileName );

  #ifdef WINDOWS
  INPUT.open( NewFileName, ios::binary|ios::nocreate );
  #else
  INPUT.open( NewFileName );
  #endif
  if( !INPUT )
       {
           GRAN[nZ1][nZ2] = NULL;
               return 0;
       }
  else
       {
           // Assign space for table
        tTABLE * pTABLE;

        pTABLE = new tTABLE;
        if( !pTABLE )
         {
          cerr << endl;
          cerr << "# Error: TBL::LoadASCII, Can't allocate " << sizeof(tTABLE)
               << " bytes !! " <<endl;
          cerr << endl;
          exit(MEMORY_ERROR);
         };
        GRAN[ nZ1 ][ nZ2 ] = pTABLE;
#ifdef GCC3
        INPUT.read( (char *) Cabecera, HEADERSIZE );
#else
        INPUT.read( (unsigned char *) Cabecera, HEADERSIZE );
#endif
    for( index_E = 0;     index_E < N_ENERGY;    index_E++    )
     for( index_PHI = 0;   index_PHI < N_PHI;     index_PHI++   )
      for( index_THETA = 0; index_THETA < N_THETA;  index_THETA++ )
       for( index_ALPHA = 0; index_ALPHA < N_ALPHA;   index_ALPHA++ )
        for( index_S = 0;     index_S     < N_S;       index_S++      )
         for( index_V = 0;     index_V     < NV3D;       index_V++      )
           {
        INPUT >> Value;
        GRAN[ nZ1 ][ nZ2 ]->TABLE
           [ index_E ][ index_PHI ][ index_THETA ][index_ALPHA]
           [ index_S ][ index_V ]     = (double ) Value;
           };
        cerr << "# Found ... " << NewFileName << endl;
    INPUT.close();
    return 1;
       }
}

///////////////////////////////////////////////////////////////////////////

void TBL::SaveBINARY( int nZ1, int nZ2, char * FileName )
{
 char NewFileName[255];
 ofstream OUTPUT;
 union { short int w; char b; } endian;
 endian.w = 0x0100;
 
 if(endian.b==1)    // BIG Endian
  {
      FileName[ 8] = 'b';
      FileName[ 9] = 'i';
      FileName[10] = 'g';   // Trying HP/UX or SUN   architecture
  }
 else
  {
      FileName[ 8] = 'l';
      FileName[ 9] = 'i';
      FileName[10] = 't';   // Trying ALPHA or INTEL architecture
  };

  TABLEDIR = getenv( "TABLES3DDIR" );  
  if( TABLEDIR != NULL ) strcpy( NewFileName, TABLEDIR );
  else                   strcpy( NewFileName, TABLEDIRDEFAULT );
  strcat( NewFileName, FileName );

  sprintf( Cabecera,"nE=%2d,n_S=%2d,nV=%1d",N_ENERGY,N_S,NV3D );

  #ifdef WINDOWS
  OUTPUT.open( NewFileName, ios::binary );
  #else
  OUTPUT.open( NewFileName );
  #endif
  if( OUTPUT )
    {
#ifdef GCC3
      OUTPUT.write( (char*) Cabecera, HEADERSIZE );
      OUTPUT.write( (char*) GRAN[nZ1][nZ2]->TABLE, sizeof(tTABLE));
#else
      OUTPUT.write( (unsigned char*) Cabecera, HEADERSIZE );
      OUTPUT.write( (unsigned char*) GRAN[nZ1][nZ2]->TABLE, sizeof(tTABLE));
#endif
      OUTPUT.close();
    }
  else 
    { 
     cerr << "# Error: TBL::SaveBinary" << endl;
     exit(WRITE_ERROR);
    };
}

///////////////////////////////////////////////////////////////////////////

void TBL::SaveASCII( int nZ1, int nZ2, char * FileName )
{
 int index_E, index_Phi, index_Theta, index_Alpha, index_S, index_V;
 double Value;
 ofstream OUTPUT;

 #ifdef WINDOWS
 OUTPUT.open( FileName, ios::binary );
 #else
 OUTPUT.open( FileName );
 #endif
 if( !OUTPUT  )
  {
   cerr << "# Error: TBL::SaveASCII" << endl;
   exit(WRITE_ERROR);
  };

 sprintf( Cabecera,"nE=%2d,n_S=%2d,nV=%1d",N_ENERGY,N_S,NV3D );
#ifdef GCC3
 OUTPUT.write( (char*) Cabecera, HEADERSIZE );
#else
 OUTPUT.write( (unsigned char*) Cabecera, HEADERSIZE );
#endif

 for( index_E = 0;   index_E < N_ENERGY;  index_E++ )
  for( index_Phi = 0; index_Phi < N_PHI;   index_Phi++ )
   for( index_Theta=0; index_Theta<N_THETA; index_Theta++ )
    for( index_Alpha = 0; index_Alpha < N_ALPHA; index_Alpha++ )
     for( index_S = 0;   index_S < N_S;     index_S++ )
     {
      for( index_V = 0;   index_V < NV3D;          index_V++ )
       {
     Value =  GRAN[ nZ1 ][ nZ2 ]->
         TABLE[ index_E ][index_Phi][index_Theta]
              [ index_Alpha ][ index_S ][ index_V ];
     OUTPUT <<  Value << " ";
       };
      OUTPUT << endl;
     };

 OUTPUT.close();
}

// Interpolate variables from Table ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void TBL::InterpolateAll_( int Atom1, int Atom2,
                 double Value_E, double Value_S,
                 double &rX1, double &rX2, double &rCOS2,
                 double &rE1, double &rE2, double &rQ,
                 double &rQElect )
{
 double V0,  V1, V16, V17;
 double e1,e2,s1,s2;
 int    Index_E, Index_S;
 int    Index;
 TableVariable Variable;

 Atom1 = Atom1 - 1;
 Atom2 = Atom2 - 1;

 if( GRAN[ Atom1 ][ Atom2 ]==NULL ) CalculeTable( Atom1+1, Atom2+1);
    
 Index_E = GetEnergyIndex( Value_E ); GetEnergyValue( Index_E, e1, e2 );
 Index_S = GetSIndex     ( Value_S ); GetSValue     ( Index_S, s1, s2 );

 if( Index_E > N_ENERGY-2 || Index_S > N_S-2 )
  {
   //MDSolver.SolveAll( );
   cerr << "# Error: Maximum energy tabulated reached\n";
   cerr << "# E= " << Value_E << " eV, p=" << Value_S << " angstroms\n";
   cerr << "# nE= " << Index_E << ", np=" << Index_S << "\n";
   exit(OUT_OF_RANGE_ERROR);
  }
 else
  {
   double K  = ( Value_E     -     e1)/( e2     -     e1 );
   double K5 = ( Value_S     -     s1)/( s2     -     s1 );

   rQ = rQElect = 0; /* */

   for ( Index = 0; Index< 5 /*NV3D*/ ; Index++ )
   {
     Variable = (TableVariable) Index;

          V0 = GRAN[ Atom1    ][ Atom2 ]->TABLE
            [ Index_E     ]
            [ 0 ][ 0 ][ 0 ]
            [ Index_S     ][ Variable ];
          V1 = GRAN[ Atom1    ][ Atom2 ]->TABLE
            [ Index_E     + 1 ]
            [ 0 ][ 0 ][ 0 ] 
            [ Index_S     ][ Variable ];
          V16= GRAN[ Atom1    ][ Atom2 ]->TABLE
            [ Index_E     ]
            [ 0 ][ 0 ][ 0 ]
            [ Index_S     + 1 ][ Variable ];
          V17= GRAN[ Atom1    ][ Atom2 ]->TABLE
            [ Index_E     + 1 ]
            [ 0 ][ 0 ][ 0 ]
            [ Index_S     + 1 ][ Variable ];

          // Interpolate in ENERGY dimension

          V0  = V0  + (V1 -V0 )*K;
          V16 = V16 + (V17-V16)*K;

          // Interpolate in S      dimension

          V0  = V0  + (V16-V0 )*K5;

    switch ( Variable )
    {
     case _X1   : rX1   = (V0>0)?0.0:V0; break;
     case _X2   : rX2   = V0;        break;
     case _COS2 : rCOS2 = V0;        break;
     case _E1   : rE1   = (V0>0)?V0:0.0;
              rE1   = (rE1<Value_E)?rE1:Value_E;
              break;
     case _E2   : rE2   = (V0>0)?V0:0.0; break;
     case _Q    : rQ    = (V0>0)?V0:0; break;
     case _QElect : rQElect = (V0>0)?V0:0; break;
    }
   }
  }
}

///////////////////////////////////////////////////////////////////////////

void TBL::CalculeTable( int A1, int A2 )
{ 
     int i, Z1, Z2;
     double W1, W2;
     AtomDefinition A;
     char Nombre[256],PNombre[256];
     char DensFile[200];
     
     A = Atoms.Search( A1 )->Data;
     Z1 = A.Z;
     W1 = A.W;
     A = Atoms.Search( A2 )->Data;
     Z2 = A.Z;
     W2 = A.W;
    
     GetSSFName(Z1,Z2,Nombre);     
     cerr << "# Generating nuclear stopping table " << Nombre << endl;

     TABLEDIR = getenv( "TABLES3DDIR" );
     if( TABLEDIR != NULL ) strcpy( PNombre, TABLEDIR );
     else                   strcpy( PNombre, TABLEDIRDEFAULT );
     strcat( PNombre, Nombre );
     sprintf( DensFile, "none" );     
     cerr << "# Z1 = " << Z1 << ",Z2 = " << Z2 
          << ",W1 = " << W1 << ",W2 = " << W2 <<endl;
     MDSolver.Init( Z1, Z2, W1, W2, DensFile, N_ENERGY, 0 );
     for( i=0; i<N_ENERGY; i++ ) MDSolver.DOIT( i );
     cerr << endl;
     MDSolver.SAVE_ASCII( PNombre, N_ENERGY );     

     LoadASCII(A1-1,A2-1,Nombre);
}
}
///////////////////////////////////////////////////////////////////////////
