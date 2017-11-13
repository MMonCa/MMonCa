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
// COL3D.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "col3d.h"

using std::endl;
using std::cerr;
using std::ofstream;
using std::cout;

namespace BCA {
// Inicializacion /////////////////////////////////////////////////////////
    
void TC::Init( double Z1_, double Z2_, double M1_,
       double M2_, const std::string &fn, int size,
       int Initialize, int SaveIt) // (YES)  NO = 0 YES = Other
{
    Out3D = 0;
    Z1 = Z1_;
    Z2 = Z2_;
    M1 = M1_*proton_mass;
    M2 = M2_*proton_mass;
    au=0.8854*bohr_radius/(pow(Z1,0.23)+pow(Z2,0.23));
    rx1 = rx2 = ry1 = ry2 = rz1 = rz2 = 0;
    Screen.Init( (int)Z1, (int)Z2, (int)SpecificPotential[(int)Z1][(int)Z2]);
    if(Initialize)
    {
      EDTable.ReadTable3D(fn, SaveIt );
    }

    // Allocate memory
    TABLE =(float*) malloc( size*N_THETA*N_PHI*N_ALPHA*N_S*NV3D*sizeof(float) );
    if( TABLE == NULL )
       {
        cerr << endl;
        cerr << "# Error: TC::Init, Can not allocate "
         << size*N_THETA*N_PHI*N_ALPHA*N_S*NV3D*sizeof(float)
         << " bytes" << endl << endl;
        exit(MEMORY_ERROR);
       };

     int ik;
     for(ik=0;ik<size*N_THETA*N_PHI*N_ALPHA*N_S*NV3D;ik++) TABLE[ik] = 1.0;
     KL  = 9e9*Z1*Z2*electron_charge*electron_charge;
 }

///////////////////////////////////////////////////////////////////////////
// Solve the trajectorie //////////////////////////////////////////////////

void TC::MakeTrajectorie( int PAINT )
{
 double rx=0,ry=0,rz=0, Rold=0, F=0;
 double fi=0,fip=0;

long numpasos=0;

Rold = 100.0e-10;
for(;;)
{
 // Trajectory solution

 rx=rx1-rx2;
 ry=ry1-ry2;
 rz=rz1-rz2;
 R=sqrt(rx*rx+ry*ry+rz*rz); // Distance between atoms
 rx/=R;
 ry/=R;
 rz/=R;             // (rx,ry,rz) Unitary vector

 // Exit condition
 if( (R>Rold)&&(R>10.0e-10) ) break;
 Rold = R;

 // vr : relative velocity between particles
 double vxx = vx1-vx2;
 double vyy = vy1-vy2;
 double vzz = vz1-vz2;

 Screen.Screening( R*1e10, fi, fip );
 fip *= 1e10;
 
 // F   : ZBL power
 F=KL/R*(fip-fi/R);

 // Inelastic losses
 // double FF3 = FFF(R,Z1,Z2);  // (S.I.) J/m
 double FF3 = 0; // No inelastico
 // Electronic stopping
 /*
  Fr=Stopping1ZBL( vr/v0, Z1, EDTable.GetRs3D(RX*1e10,RY*1e10,RZ*1e10) );
  double FF = Fr*electron_charge/ sqrt(vx1*vx1+vy1*vy1+vz1*vz1);
 */
 double FF = 0;

 // Accelerations
 ax1=(-F*rx-FF*vx1 -FF3*vxx )/ M1;
 ay1=(-F*ry-FF*vy1 -FF3*vyy )/ M1;
 az1=(-F*rz-FF*vz1 -FF3*vzz )/ M1;

 double FE = F/M2;
 ax2= FE*rx;
 ay2= FE*ry;
 az2= FE*rz;

 // New positions
 double r1 = vx1*h+ax1*h2h;
 double r2 = vy1*h+ay1*h2h;
 double r3 = vz1*h+az1*h2h;

 rx1+=r1;
 ry1+=r2;
 rz1+=r3;

 rx2+=vx2*h+ax2*h2h;
 ry2+=vy2*h+ay2*h2h;
 rz2+=vz2*h+az2*h2h;

 // New velocities
 vx1+=hh*(ax1+ax10);
 vy1+=hh*(ay1+ay10);
 vz1+=hh*(az1+az10);

 vx2+=hh*(ax2+ax20);
 vy2+=hh*(ay2+ay20);
 vz2+=hh*(az2+az20);

 // Add inelastic losses
 //q+=FF3*sqrt(vxx*vxx+vyy*vyy+vzz*vzz)*sqrt(r1*r1+r2*r2+r3*r3); // J

 // Add electronic losses
 //qElect += Fr*sqrt(r1*r1+r2*r2+r3*r3);

 // copia de aceleraciones
 ax10=ax1; ay10=ay1; az10=az1;
 ax20=ax2; ay20=ay2; az20=az2;

 numpasos++;

 }; // End of FOR loop

 q = q / electron_charge; // Put in correct units (eV)
 qElect = qElect * 1e10;
}

///////////////////////////////////////////////////////////////////////////

void TC::SetParamPSystem( double ENERGY, double Xx, double Yy, double Zz )
{
     // Set initial positions in laboratory system
      // Do the transformation of coordinates. From ('') to (w ')
      double X = MI[0][0] * Xx + MI[0][1] * Yy + MI[0][2] * Zz;
      double Y = MI[1][0] * Xx + MI[1][1] * Yy + MI[1][2] * Zz;
      double Z = MI[2][0] * Xx + MI[2][1] * Yy + MI[2][2] * Zz;
      // Projectile
      rx1 = X;  // In meters
      ry1 = Y;
      rz1 = Z;
      // Target. No transformation needed.
      rx2 = 0.;
      ry2 = 0.;
      rz2 = 0.;

      h=0.001e-10/sqrt(2.*ENERGY/M1);  // Time in seconds
      hh  = h * 0.5;
      h2h = h*h * 0.5;

      // Velocity
      double V = sqrt(2.*ENERGY/M1);

      q =0.;    //    ap=1.;
      qElect = 0;

     // Set initial velocities
      // Projectile
      // Do the transformation of coordinates. From ('') to (w ')
      vx1 = MI[0][0] * V; // + MI[0][1] * 0 + MI[0][2] * 0;
      vy1 = MI[1][0] * V; // + MI[1][1] * 0 + MI[1][2] * 0;
      vz1 = MI[2][0] * V; // + MI[2][1] * 0 + MI[2][2] * 0;
      // Target. No transformation needed.
      vx2 = 0.;
      vy2 = 0.;
      vz2 = 0.;
      // Clear accelerations
      ax1 =0.; ay1 =0.; az1 =0.;
      ax10=0.; ay10=0.; az10=0.;
      ax2 =0.; ay2 =0.; az2 =0.;
      ax20=0.; ay20=0.; az20=0.;
}

///////////////////////////////////////////////////////////////////////////

void TC::MakeMatrix( double sinTHETA, double cosTHETA,
             double sinPHI, double cosPHI,
               double sinALPHA, double cosALPHA )
{
      // The inverse matrix is equal to the trasposed matrix
      // Make the transformation matrix

      MI[0][0] = MD[0][0] =   cosTHETA*cosPHI*cosALPHA - sinPHI*sinALPHA;
      MI[0][1] = MD[1][0] = - sinTHETA*cosALPHA;
      MI[0][2] = MD[2][0] =   cosTHETA*sinPHI*cosALPHA + cosPHI*sinALPHA;

      MI[1][0] = MD[0][1] =   sinTHETA*cosPHI;
      MI[1][1] = MD[1][1] =   cosTHETA;
      MI[1][2] = MD[2][1] =   sinTHETA*sinPHI;

      MI[2][0] = MD[0][2] = - cosTHETA*cosPHI*sinALPHA - sinPHI*cosALPHA;
      MI[2][1] = MD[1][2] =   sinTHETA*sinALPHA;
      MI[2][2] = MD[2][2] = - cosTHETA*sinPHI*sinALPHA + cosPHI*cosALPHA;

} // End of MakeMatrix

///////////////////////////////////////////////////////////////////////////

double TC::GetResultPSystem( double &x1, double &x2, double &cos1, double &cos2,
             double &E1, double &E2, double &Q, double &QElect )
{
     double tan2,v1;

     // Solve the ENERGIES
     E1=0.5*M1*(vx1*vx1+vy1*vy1+vz1*vz1);
     E2=0.5*M2*(vx2*vx2+vy2*vy2+vz2*vz2);

     // Make the transformation
     double Vx1 = MD[0][0] * vx1 + MD[0][1] * vy1 + MD[0][2] * vz1;
     double Vy1 = MD[1][0] * vx1 + MD[1][1] * vy1 + MD[1][2] * vz1;
     double Vz1 = MD[2][0] * vx1 + MD[2][1] * vy1 + MD[2][2] * vz1;

     double X2  = MD[0][0] * rx2 + MD[0][1] * ry2 + MD[0][2] * rz2;
     double Z2  = MD[2][0] * rx2 + MD[2][1] * ry2 + MD[2][2] * rz2;

     double Vx2 = MD[0][0] * vx2 + MD[0][1] * vy2 + MD[0][2] * vz2;
     double Vy2 = MD[1][0] * vx2 + MD[1][1] * vy2 + MD[1][2] * vz2;
     double Vz2 = MD[2][0] * vx2 + MD[2][1] * vy2 + MD[2][2] * vz2;

double VV = sqrt( vx1*vx1+vy1*vy1+vz1*vz1 );
double c1 = vx1/VV*MI[0][0] + vy1/VV*MI[1][0] + vz1/VV*MI[2][0];
double a  = rx1*MI[0][0] + ry1*MI[1][0] + rz1*MI[2][0];
double b  = sqrt(  rx1*rx1+ry1*ry1+rz1*rz1 - a*a );
double ttan = tan( acos(c1) );
if( fabs(ttan) > 1e-10 ) x1 =  (a - (b-S_3D)/ttan);
else             x1 = 0.0;          // 19991111

     if( fabs( Vx2 )<1e-10 )
       if( Vz2>0 ) tan2 =  1e30; else tan2 = -1e30;
     else tan2 = Vz2/Vx2;
     if(fabs(tan2)>1e-10) x2=X2-Z2/tan2;
     else             x2=X2;

     v1=sqrt(Vx1*Vx1+Vy1*Vy1+Vz1*Vz1); if (v1<1e-5) v1=1e-5;
     cos1=Vx1/v1;

     v1=sqrt(Vx2*Vx2+Vy2*Vy2+Vz2*Vz2); if (v1<1e-5) v1=1e-5;
     cos2=Vx2/v1;

     x1 = x1 * 1e10;
     x2 = x2 * 1e10;
     E1 = E1 / electron_charge;
     E2 = E2 / electron_charge;
     Q  = q;
     QElect = qElect;
     return q;
}

///////////////////////////////////////////////////////////////////////////

void TC::SolveAll(
    double ENERGY, double S, double Ghi, double THETA, double PHI, double ALPHA,
    double &x1, double &x2, double &cos2, double &E1, double &E2, double &Q,
    double &QElect )
{
  double Xx, Yy, Zz, cos1;

  ENERGY = ENERGY * electron_charge;
  S  = S   * 1e-10; // In metres
  Ghi    = Ghi * 1e-10; // In metres
  MakeMatrix(sin(THETA),cos(THETA), sin(PHI), cos(PHI), sin(ALPHA), cos(ALPHA));
  S_3D = S;
  Xx = -3.0 ; // - Ghi;
  Yy = 0.0;
  Zz = S;
  SetParamPSystem( ENERGY, Xx, Yy, Zz );
  MakeTrajectorie();
  GetResultPSystem( x1,x2,cos1,cos2,E1,E2,Q, QElect );
}

///////////////////////////////////////////////////////////////////////////

void TC::MAKEOUT( int i1, int i2, int i3, int i4, int i5, int Print)
{
     double E1,E2,cos1,cos2,x1,x2, Q, QElect;
     float * pTABLE;

     GetResultPSystem( x1, x2, cos1, cos2, E1, E2, Q, QElect );
     pTABLE = TABLE + i1 * NV3D * N_S * N_ALPHA * N_THETA * N_PHI
            + i2 * NV3D * N_S * N_ALPHA * N_THETA
            + i3 * NV3D * N_S * N_ALPHA
            + i4 * NV3D * N_S
            + i5 * NV3D;
     *(pTABLE     ) = x1 ;
     *(pTABLE + 1 ) = x2 ;
     *(pTABLE + 2 ) = cos2;
     *(pTABLE + 3 ) = E1 ;
     *(pTABLE + 4 ) = E2 ;
     *(pTABLE + 5 ) = Q;
     *(pTABLE + 6 ) = QElect;    
     
    if (Print)
     {
       cout << endl;
       cout << x1   << " " << x2    << " "  << cos1 << " " << cos2  << " "
        << E1   << " " << E2    << " "  << Q    << " "
        << QElect     << endl << endl;
     };
}

///////////////////////////////////////////////////////////////////////////

void TC::SAVE_ASCII( char* pfnt, int size )
{
      char Cabecera[HEADERSIZE+1];
      ofstream OUT;
      int i1, i2, i3, i4, i5;
      float * pT;
      
      for(i1=0;i1<HEADERSIZE-1;i1++)Cabecera[i1]='%';
      Cabecera[HEADERSIZE-1]='\n';
      sprintf(Cabecera, "nE=%2d,nS=%2d,nV=%1d",N_ENERGY,N_S,NV3D );

      #ifdef WINDOWS
      OUT.open( pfnt, ios::binary );
      #else
      OUT.open( pfnt );
      #endif
      if( !OUT )
       { cerr << "# Error: TC::SAVE_ASCII, Can't open file " 
              << pfnt << endl; 
         exit(READ_ERROR);
       };
      
      OUT.setf(std::ios::fixed);
      OUT.precision(6);
#ifdef GCC3
      OUT.write( (char*) &Cabecera, HEADERSIZE );
#else
      OUT.write( (unsigned char*) &Cabecera, HEADERSIZE );
#endif
      pT = TABLE;
      for( i1=0; i1<size;     i1++ )
      for( i2=0; i2<N_PHI;    i2++ )
      for( i3=0; i3<N_THETA;  i3++ )
      for( i4=0; i4<N_ALPHA;  i4++ )
      for( i5=0; i5<N_S;      i5++ )
         {
          pT = TABLE
            + i1 * NV3D * N_S * N_ALPHA * N_THETA * N_PHI
            + i2 * NV3D * N_S * N_ALPHA * N_THETA
            + i3 * NV3D * N_S * N_ALPHA
            + i4 * NV3D * N_S
            + i5 * NV3D ;
                
      OUT.width( 15 );
      OUT <<  *(pT     ) << " ";
      OUT.width( 15 );
      OUT <<  *(pT + 1 ) << " ";
      OUT.width( 15 );
      OUT <<  *(pT + 2 ) << " ";
      OUT.width( 15 );  
      OUT <<  *(pT + 3 ) << " ";
      OUT.width( 15 );
      OUT <<  *(pT + 4 ) << " ";
      OUT.width( 15 );
      OUT <<  *(pT + 5 ) << " ";
      OUT.width( 15 );
      OUT <<  *(pT + 6 ) << endl;

          };
      OUT.close();
}

///////////////////////////////////////////////////////////////////////////

void TC::SAVE_BIN( char *pfnt, int size )
{
     char Cabecera[HEADERSIZE+1];
     ofstream FF;

     sprintf(Cabecera, "nE=%2d,nS=%2d,nV=%1d",N_ENERGY,N_S,NV3D );

     #ifdef WINDOWS
     FF.open(pfnt,ios::binary);
     #else
     FF.open(pfnt);
     #endif
#ifdef GCC3
     FF.write( (char*) &Cabecera, HEADERSIZE );
     FF.write( (char*) TABLE,
      size * N_S * N_ALPHA * N_THETA * N_PHI * NV3D * sizeof( float ) );
#else     
     FF.write( (unsigned char*) &Cabecera, HEADERSIZE );
     FF.write( (unsigned char*) TABLE,
      size * N_S * N_ALPHA * N_THETA * N_PHI * NV3D * sizeof( float ) );
#endif     
     FF.close();
}

///////////////////////////////////////////////////////////////////
// Not used
void TC::DOIT_PVM( int i_ENERGY, int i_PHI, int i_E, int i_P )
    {
      //
      // Define all the loops for the input parameters to do the whole range
      //
      double Xx, Yy, Zz;
      double S, ALPHA, PHI, THETA, ENERGY, ENERGY_, Temp;
      double cosPHI, sinPHI, cosTHETA, sinTHETA;

      time_t T0, T2;
      T0 = time( NULL );
      GetEnergyValue( i_ENERGY, ENERGY_, Temp );
      ENERGY = ENERGY_ * electron_charge; // In Coulomb
      GetPhiValue( i_PHI, PHI, Temp );
      cosPHI = cos(PHI);
      sinPHI = sin(PHI);
      int i_THETA = 0;
      //for(int i_THETA=0; i_THETA < N_THETA; i_THETA++ )
      //{
        GetThetaValue( i_THETA, THETA, Temp );
        sinTHETA = sin(THETA);
        cosTHETA = cos(THETA);
        int i_ALPHA = 0;
        //for(int i_ALPHA=0; i_ALPHA < N_ALPHA; i_ALPHA++ )
        //{
        	GetAlphaValue( i_ALPHA, ALPHA, Temp );
        	MakeMatrix( sinTHETA, cosTHETA,
                sinPHI,   cosPHI,
                sin(ALPHA), cos(ALPHA) );

			for(char i_S = 0; i_S < N_S; i_S++) //no idea why int causes wrong loop optmization
			{
			   GetSValue( i_S, S, Temp );
			   S = S * 1e-10; // In metres
			   S_3D = S;
			   Yy = 0.0;  // In meters
			   Zz = S; // In meters
			   Xx  = - 3.0e-10;  // In meters
			   SetParamPSystem( ENERGY, Xx, Yy, Zz );
			   MakeTrajectorie();
			   MAKEOUT( i_E, i_P, i_THETA, i_ALPHA, i_S );
           }
       //}
     //}
     T2 = time( NULL );
     char MSG[80];
     sprintf( MSG,"(%3d:%3d) %lf seconds", i_ENERGY,i_PHI,
         difftime( T2, T0 ) );
         cerr << MSG;
         for(unsigned i=0; i<strlen(MSG); i++ ) cerr << "\b";
 } // End of DOIT()

///////////////////////////////////////////////////////////////////////////

void TC::DOIT( int i_ENERGY )
{
  char MSG[81];
  int i_PHI;
  time_t T0, T2;
  T0 = time( NULL );
  for( i_PHI=0;   i_PHI< N_PHI;     i_PHI++   )
   {
     DOIT_PVM( i_ENERGY, i_PHI, i_ENERGY, i_PHI );
   };
  T2 = time( NULL );
  sprintf( MSG,"# Index %3d/%3d calculated in %.2lf seconds", i_ENERGY,
        N_ENERGY-1, difftime( T2, T0 ) );
  cerr << MSG;
  for(unsigned i=0; i<strlen(MSG); i++ ) cerr << "\b";
} // End of DOIT()


///////////////////////////////////////////////////////////////////

void TC::DOITTEST( int i_ENERGY, int i_PHI, int i_THETA,
           int i_ALPHA, int i_S, int Print )
{
      double Xx, Yy, Zz;
      double S, ALPHA, PHI, THETA, ENERGY, Temp;

      GetEnergyValue( i_ENERGY, ENERGY, Temp );
      GetPhiValue   ( i_PHI   , PHI,    Temp );
      GetThetaValue ( i_THETA , THETA,  Temp );
      GetAlphaValue ( i_ALPHA,  ALPHA,  Temp );
      GetSValue ( i_S     , S,      Temp );

      ENERGY = ENERGY * electron_charge;
      S  = S   * 1e-10;
      MakeMatrix( sin(THETA), cos(THETA),
              sin(PHI),   cos(PHI),
              sin(ALPHA), cos(ALPHA) );
      S_3D = S;
      Xx   = - 3e-10;    // In meters
      Yy   = 0.0;        // In meters
      Zz   = S;      // In meters

      // Make the calculation

     #ifdef TEST
      cout <<    " Ghi = " << Ghi * 1e10;
      cout << " (A), S = " << S   * 1e10;
      cout << " (A), E = " << ENERGY / e << " eV" << endl;
     #endif
      SetParamPSystem( ENERGY, Xx, Yy, Zz );
      MakeTrajectorie();
      MAKEOUT( i_ENERGY, i_PHI, i_THETA, i_ALPHA, i_S, Print );
    }

static TC MDSolver;
}
