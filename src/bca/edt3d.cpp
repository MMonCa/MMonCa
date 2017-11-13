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
// EDT 3D.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "edt3d.h"

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;

namespace BCA {

int EDT::ReadTable3D(const std::string &FileName, int SaveIt = 1 )
{
 int I, J, K;
 double Value;
 char ASCII[80], NewFileName[255];
 ifstream FP;
 ofstream OF;

 union { short int w; char b; } endian; // endian.b==1 => Big    Endian
 endian.w = 0x0100;                     // endian.b==0 => Little Endian

 char *EDTDIR, NewFileName2[255];
 EDTDIR = getenv( "EDTDIR" ); // Read the environment variable
 if( EDTDIR != NULL ) strcpy( NewFileName2, EDTDIR );
 else                 strcpy( NewFileName2, EDTDIRDEFAULT );
 strcat( NewFileName2, FileName.c_str());
                                        
 // Try to open binary file
 strcpy( NewFileName, NewFileName2 );
 if( endian.b ) strcat( NewFileName, ".big" );  // BIG  Endian arch.
 else           strcat( NewFileName, ".lit" );  // LITTLE Endian arch.
#ifdef WINDOWS
 FP.open( NewFileName, ios::in|ios::binary|ios::nocreate );
#else
 FP.open( NewFileName );
#endif

 if( FP )
  {
      char car;
      
      FP.read( &car,1 ); while(car!=0x0a) { FP.read( &car,1 ); }
      FP.read( &car,1 ); while(car!=0x0a) { FP.read( &car,1 ); }
      FP >> LADOX >> N_POINTSX;
      FP >> LADOY >> N_POINTSY;      
      FP >> LADOZ >> N_POINTSZ;
      FP.read( &car,1 ); while(car!=0x0a) { FP.read( &car,1 ); }

      SizeOfDATA = N_POINTSX*N_POINTSY*N_POINTSZ*sizeof(double);
      DATA = (double*) malloc( SizeOfDATA );
      if( DATA==NULL ) { perror("# No memory"); exit(MEMORY_ERROR); };
#ifdef GCC3
      FP.read( (char*) DATA, SizeOfDATA );
#else      
      FP.read( (unsigned char*) DATA, SizeOfDATA );
#endif
      FP.close();
      cerr << "# 3D binary density table read\n#\t(" << NewFileName << ") "
           << FP.gcount() << "/" << SizeOfDATA << endl;
  }
 else // Else try to open ASCII file
  {
  #ifdef WINDOWS
   FP.open(NewFileName2,ios::nocreate);
  #else
   FP.open(NewFileName2);
  #endif 
    if (!FP)
     { // Create 3D file
      cerr << "# Error: EDT::ReadTable3D, File not found: " << NewFileName2 << endl;
      return 1;
     };
        char car;
        
        FP.read( &car,1 ); while(car!=0x0a) { FP.read( &car,1 ); }  
        FP.read( &car,1 ); while(car!=0x0a) { FP.read( &car,1 ); }  
        FP >> LADOX >> N_POINTSX;
        FP >> LADOY >> N_POINTSY;      
        FP >> LADOZ >> N_POINTSZ;
        FP.read( &car,1 ); while(car!=0x0a) { FP.read( &car,1 ); }  

        SizeOfDATA = N_POINTSX*N_POINTSY*N_POINTSZ*sizeof(double);
        DATA = (double*) malloc( SizeOfDATA );
        if( DATA==NULL ) { perror("# No memory"); exit(MEMORY_ERROR); };

    double* pDATA = DATA; 
     
    for( I=0; I<N_POINTSX; I++ )
     for( J=0; J<N_POINTSY; J++ )
      for( K=0; K<N_POINTSZ; K++ )
       {
         FP >> ASCII;
         Value = atof( ASCII );
         *(pDATA++ ) = Value;
       };
    FP.close();
   cerr << "# 3D ASCII density table read\n\t(" << NewFileName2 << ")" 
    << endl;
   if( SaveIt) SaveBinTable3D( NewFileName2 );
  };

 dX = N_POINTSY*N_POINTSZ;
 dY = N_POINTSZ;

 return 0; 
}

////////////////

void EDT::InitEmptyTable( double LX, double LY, double LZ,
              int nX, int nY, int nZ )
{
 N_POINTSX = nX;
 N_POINTSY = nY;
 N_POINTSZ = nZ;
 LADOX = LX;
 LADOY = LY;
 LADOZ = LZ;

 dX = nY*nZ;
 dY = nZ;
 SizeOfDATA = nX*nY*nZ*sizeof(double);
 DATA = (double*) malloc( SizeOfDATA );
 if( DATA==NULL ) { perror("# No memory"); exit(MEMORY_ERROR); };

}

///////////////

void EDT::SetValue( int i, int j, int k, double Value )
{
 *(DATA + i*dX + j*dY + k) = Value;
}

///////////////////////////////////////////////////////////////////////////

inline double EDT::Interpolate3D( double X, double Y, double Z )
{ // Coordinates are in angstroms
   int x0, y0, z0;
   double X0, Y0, Z0;
   double Value;
   double V1, V2, V3, V4, V5, V6, V7, V8;
   double V12, V34, V56, V78;
   double V14, V58;
  //
  // The density table has a LADO periodicity
  //

   X = X + LADOX/2.;
   if ( X < 0 ) X = LADOX - LADOX*( - X/LADOX + (int) (X/LADOX) );
   else         X =         LADOX*(   X/LADOX - (int) (X/LADOX) );
   X = X - LADOX/2.;

   Y = Y + LADOY/2.;
   if( Y < 0 )  Y = LADOY - LADOY*( -Y/LADOY + (int) (Y/LADOY) );
   else         Y =         LADOY*(  Y/LADOY - (int) (Y/LADOY) );
   Y = Y - LADOY/2.;

   Z = Z + LADOZ/2.;
   if( Z < 0 ) Z = LADOZ - LADOZ*( -Z/LADOZ + (int) (Z/LADOZ) );
   else        Z =         LADOZ*(  Z/LADOZ - (int) (Z/LADOZ) );
   Z = Z - LADOZ/2.;

   X = X *N_POINTSX/LADOX;
   Y = Y *N_POINTSY/LADOY;
   Z = Z *N_POINTSZ/LADOZ;

   // Get the lower indexes
   x0 = (int) (X + N_POINTSX/2);
   y0 = (int) (Y + N_POINTSY/2);
   z0 = (int) (Z + N_POINTSZ/2);

   // Extrapolation
  #ifndef NOEXTRAPOLATION
   if( x0 <  0 )      x0 = 0;
   if( x0 >= N_POINTSX-1 ) x0 = N_POINTSX-2;

   if( y0 <  0 )      y0 = 0;
   if( y0 >= N_POINTSY-1 ) y0 = N_POINTSY-2;

   if( z0 <  0 )      z0 = 0;
   if( z0 >= N_POINTSZ-1 ) z0 = N_POINTSZ-2;

   // No extrapolation
  #else
   if( (x0 <   0) || (x0 >= N_POINTSX-1) ||
       (y0 <   0) || (y0 >= N_POINTSY-1) ||
       (z0 <   0) || (z0 >= N_POINTSZ-1) )
     return (0);
  #endif

   X0 = ( x0 - N_POINTSX/2. );
   Y0 = ( y0 - N_POINTSY/2. );
   Z0 = ( z0 - N_POINTSZ/2. );

   double K  = ( X - X0 );
   double K2 = ( Y - Y0 );
   double K3 = ( Z - Z0 );

   V1 = *( DATA + (x0  )*dX + (y0  )*dY + (z0  ) );
   V2 = *( DATA + (x0+1)*dX + (y0  )*dY + (z0  ) );
   V3 = *( DATA + (x0  )*dX + (y0+1)*dY + (z0  ) );
   V4 = *( DATA + (x0+1)*dX + (y0+1)*dY + (z0  ) );
   V5 = *( DATA + (x0  )*dX + (y0  )*dY + (z0+1) );
   V6 = *( DATA + (x0+1)*dX + (y0  )*dY + (z0+1) );
   V7 = *( DATA + (x0  )*dX + (y0+1)*dY + (z0+1) );
   V8 = *( DATA + (x0+1)*dX + (y0+1)*dY + (z0+1) );

   double max = V1;
   if(V2>max) max=V2;
   if(V3>max) max=V3;
   if(V4>max) max=V4;
   if(V5>max) max=V5;
   if(V6>max) max=V6;
   if(V7>max) max=V7;
   if(V8>max) max=V8;

   double VM = 1e-10;

   if (V2<=VM) V2 = max;
   if (V1<=VM) V1 = max;
   if (V4<=VM) V4 = max;
   if (V3<=VM) V3 = max;
   if (V6<=VM) V6 = max;
   if (V5<=VM) V5 = max;
   if (V8<=VM) V8 = max;
   if (V7<=VM) V7 = max;

   // Interpolate in X dimension

   V12 = V1 + (V2-V1)*K;
   V34 = V3 + (V4-V3)*K;
   V56 = V5 + (V6-V5)*K;
   V78 = V7 + (V8-V7)*K;

   // Interpolate in Y dimension
   V14 = V12 + (V34-V12)*K2;
   V58 = V56 + (V78-V56)*K2;

   // Interpolate in Z dimension
   Value = V14 + (V58-V14)*K3;
   if(Value<=0.0) Value = VM;
   return(Value);
}

///////////////////////////////////////////////////////////////////////////

inline double EDT::NoInterpolate3D( double X, double Y, double Z )
{ // Coordinates are in angstroms
   int x0, y0, z0;
   double Value;

   // Get the lower indexes
   x0 = (int) (X*N_POINTSX/LADOX + N_POINTSX/2.);
   y0 = (int) (Y*N_POINTSY/LADOY + N_POINTSY/2.);
   z0 = (int) (Z*N_POINTSZ/LADOZ + N_POINTSZ/2.);

   // No extrapolation
   if( (x0 <   0) || (x0 >= N_POINTSX-1) ||
       (y0 <   0) || (y0 >= N_POINTSY-1) ||
       (z0 <   0) || (z0 >= N_POINTSZ-1) )
     return (0.0);

   Value = *( DATA + (x0  )*dX + (y0  )*dY + (z0  ) );
   return( Value );
}

///////////////////////////////////////////////////////////////////////////

inline double EDT::GetRs3D( double X, double Y, double Z )
{  // Coordinates are in angstroms
  double Value = Interpolate3D(X,Y,Z);
  assert( Value != 0.0 );

  return( pow( 0.75/M_PI/Value,1./3.)/bohr_radius*1e-10 );
}

///////////////////////////////////////////////////////////////////////////

void EDT::Cut3D( double x0, double y0, double z0,   // Point
         double A,  double B,  double C  )  // Direction
{
  int I, J, K;
  double XValue, YValue, ZValue, Value;
  double KKX = LADOX/N_POINTSX; // Width(A) / N_POINTS
  double KKY = LADOY/N_POINTSY;
  double KKZ = LADOZ/N_POINTSZ;
  double D = - A * x0 - B * y0 - C * z0;

  for( I=0; I<N_POINTSX; I++ )  // Scan X coordinates
   for( J=0; J<N_POINTSY; J++ ) // Scan Y coordinates
    for( K=0; K<N_POINTSZ; K++ ) // Scan Z coordinates
     {
      // The double coordinates in (A)
      XValue = (I-N_POINTSX/2) * KKX;
      YValue = (J-N_POINTSY/2) * KKY;
      ZValue = (K-N_POINTSZ/2) * KKZ;

      // Evaluate the plane equation
      Value = A * XValue + B * YValue + C * ZValue + D;

      // The condition
      if( Value > 0.0 ) 
        *( DATA + (I  )*dX + (J  )*dY + (K  ) ) = 1e-30;
     };
}

///////////////////////////////////////////////////////////////////////////

void EDT::SaveASCIITable3D( char * FileName )
{
 int I, J, K;
 double Value;
 ofstream FP;

 FP.open(FileName);
 if (!FP)
  {
   cerr << "# Error: EDT::SaveASCIITable3D, File not found: " << FileName << endl;
   exit(WRITE_ERROR);
  };
  
 cerr << "# ASCII 3D density table"; 

 FP << " EDT File (" << FileName << ") for ?" << endl;
 FP << "  LADO (A)   N_POINTS" << endl;
 FP << LADOX << " " << N_POINTSX << endl;
 FP << LADOY << " " << N_POINTSY << endl;
 FP << LADOZ << " " << N_POINTSZ << endl;
  
 double* pDATA = DATA;

 for( I=0; I<N_POINTSX; I++ )
  for( J=0; J<N_POINTSY; J++ )
   for( K=0; K<N_POINTSZ; K++ )
     {
      Value = *(pDATA++);
      FP << Value << endl;
     };
 FP.close();
 cerr << " written (" << FileName << ")" << endl;
}

// Save the table in architecture dependent binary form
void EDT::SaveBinTable3D( char * FileName )
{
 ofstream FP;
 char NewFileName[255];

 union { short int w; char b; } endian; // endian.b==1 => Big    Endian
 endian.w = 0x0100;                     // endian.b==0 => Little Endian

 strcpy( NewFileName, FileName );

 if( endian.b ) strcat( NewFileName, ".big" ); // BIG    Endian arch.
 else           strcat( NewFileName, ".lit" ); // LITTLE Endian arch.
 #ifdef WINDOWS
 FP.open( NewFileName, ios::binary );
#else
 FP.open( NewFileName );
#endif
  if( FP )
    {
      char HeaderFile[301];

      for( unsigned i=0; i<sizeof(HeaderFile); i++ ) HeaderFile[i] = ' ';
      cerr << "# Binary 3D density table";
      sprintf( HeaderFile, 
         " EDT File ( %s ) \n  LADO (A)   N_POINTS\n %f %d\n %f %d\n %f %d ",
         NewFileName, LADOX, N_POINTSX, LADOY, N_POINTSY, LADOZ, N_POINTSZ );
      HeaderFile[300]=0x0a;
  
#ifdef GCC3
      FP.write( (char*) HeaderFile, sizeof(HeaderFile));
      FP.write( (char*) DATA, SizeOfDATA );
#else     
      FP.write( (unsigned char*) HeaderFile, sizeof(HeaderFile));
      FP.write( (unsigned char*) DATA, SizeOfDATA );
#endif     
      FP.close();      
    }
  else
    {
      cerr << "# Error: EDT::SaveBinTable3D, File not found: " << FileName << endl;
      exit(WRITE_ERROR);
    };
    
 cerr << " written (" << NewFileName << ")" << endl;
}

///////////////////////////////////////////////////////////////////////////
//
//  DENSITY
//

double DensityClass::LossesNew( double Z1, double M1, Vector PAbs, Vector Dir,
           double ENERGY, double E1, double AllowedDistance, int Z2 )
  {
   // Variables -----------------
   Vector DR, Delta;
   double V, Fr, Distance, q, DELTA;
   double Et;
   double au = 0.4683766/(pow((double )Z1,0.23)+pow((double )Z2,0.23));  // (A)
   double K = Z1 * Z2 * electron_charge * 9e19;
   double Epotencial, VV;
   SSF Ssf;

   // Code -------------------- 
   DR      = PAbs;                        // (A)
   AllowedDistance = AllowedDistance;     // (A)
   DELTA = AllowedDistance/NPASOS_S;
   Delta   = Dir*DELTA;                   // (A)
   Distance = 0;
   q = 0.0;
   Et= 0.0;
   double CTE1 = sqrt( CTE/M1 );
   double CTE2 = (E1-ENERGY)/AllowedDistance;
   while( Distance < AllowedDistance )
   { 
    // Potential energy calculation ---------------------
    if(Z2!=0)
    {
     Epotencial=0.0;
     LinkObject<Vector> * O;
     O=LP.PointerToFirstData();
     do
     {
      double DD = (double ) ( O->Data - DR );
      Ssf.ScreeningOnly( DD/au, VV );
      Epotencial += VV / DD;
     }
     while( O=LP.PointerToNextData(),O!=NULL);
     Epotencial = K * Epotencial; // (eV)
    }
    else Epotencial=0.0;
     
    // Electronic stopping calculus  --------------------
    Et += CTE2;
    double EEE = ENERGY - q - Epotencial -Et;
    if( EEE < 0 ) Fr = 0; 
    else
     {
       double rs = GetRs3D( DR.X,DR.Y,DR.Z);
       V  = CTE1*sqrt( EEE );               // Atomic units
       Fr = STPO.Stopping(V,Z1,rs);         // Optimized
      }
    DR = DR + Delta;
    q        += Fr;
    Distance += DELTA;
   };

   return q*DELTA;
  }
  
}

