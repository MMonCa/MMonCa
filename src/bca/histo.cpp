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
// HISTO.H
//
///////////////////////////////////////////////////////////////////////////
#include "histo.h"

using std::cerr;
using std::endl;
using std::ofstream;

namespace BCA {

// ----------------------------------------------------------

void Hist::Initialize(double MaxDepth ,double D )
   {
     int i;
     
     Dose=D;
     DMax = 1.0;
     DMin = 0.0;
        
     // Isospaced splitting depths        
     for(i=0; i<N_BIN+1; i++ ) { H[i] = 0.0; HR[i] = 0.0; };

     LVI.DeleteAll();
   }    
   
void Hist::Add( double D, double W, double O )
   {
     VirtualIon VI;
     VI.Depth  = D;
     VI.Weigth = W;
     VI.Other  = O;
     LVI.Add(VI);
   }
   
void Hist::Do()         // Solve the histogram
   {
     LinkObject<VirtualIon> * LOVI;
     double D,W,O;
     int i,iO;

     for(i=0; i<N_BIN+1; i++ ) 
      {
       H[i] = 0.0;
       HR[i] = 0.0;
       H_err[i] = 0.0;
      };    
     Outside = 0;
      
     // Find the maximum
     DMax = 0.0;
     DMin = 1e10;
     DMaxO =0.0;
     LOVI=LVI.PointerToFirstData();
     do
      {
       D = LOVI->Data.Depth;
       O = LOVI->Data.Other;
       if( D>DMax ) DMax = D;
       if( D<DMin ) DMin = D;
       if( O>DMaxO) DMaxO= O; 
       LOVI=LVI.PointerToNextData();
      }
     while(LOVI != NULL);
     // Do the histogram    
     RealIons = 0.0;
     LOVI=LVI.PointerToFirstData();
     do
      {
       D = LOVI->Data.Depth;
       O = LOVI->Data.Other;
       W = LOVI->Data.Weigth;
//       i = (int) (D/DMax*N_BIN);
       i = floor(D/DMax*N_BIN); // To prevent negative coordinate atoms to be
                                // included in bin number 0.

       iO= (int) (O/DMaxO*N_BIN);

       if(i>=0)
       {
        H[i] += W;
        HR[iO] +=W;
        H_err[i] +=1;
       }
       else Outside+=W;         // Counts outside atoms
       
       RealIons += W;
      
       LOVI=LVI.PointerToNextData();
      }
     while(LOVI != NULL);

    for(i=0;i<N_BIN;i++) 
     if(H_err[i]!=0) H_err[i] = H[i] / sqrt( H_err[i] );
     else            H_err[i] = 1e10;
   }

void Hist::SaveHisto( int NProj,
                    int HST, char* HSTFile,
                    int FullDoseOutput, int KK, unsigned long TotalRealIons )
   {
     int i;
     double K;
     ofstream OF;
     char F[300];
     
     if( FullDoseOutput ) K = Dose*N_BIN/RealIons/1e-8/DMax;
     else                 K = Dose*N_BIN/RealIons/1e-8/DMax
                                *(double )(KK)/TotalRealIons;
     
     sprintf(F,"%s/.%s_%02d",DIRECTORY,HSTFile,NProj); // Hidden file (Unix)
     OF.open( F );
     if(!OF)
        {cerr << "# Error: SaveHisto, Can't create "<<F<<"\n";exit(WRITE_ERROR);}
     OF << "# " << KK << " " << RealIons << " " << LVI.NumberOfElements() << " "
                << Outside*K << endl;
     for(i=0; i<N_BIN; i++ )
     {  
      OF << i*DMax/N_BIN*0.1 << " " << H[i]*K << " "
         << (H[i]-H_err[i])*K << " " << (H[i]+H_err[i])*K << endl;
     }  
     
// If you activate these lines the JoinDoseSplitting will not work
//     OF << "# The two first numbers after the # means:\n";
//     OF << "#  Trajectories, double ions = NumberOfImplants*EffectiveWeigth, Virtual ions\n";
//     OF << "# Next numbers are:\n";
//     OF << "#  Depth(nm), Conc. (at/cm^3), Lower error bar, Upper error bar\n";
     OF.close(); 
   }

void Hist::SaveHistoR(int NProj,int HST,char* HSTFile)
   {
     int i;
     ofstream OF;   
     char F[300];
     double K = Dose*N_BIN/RealIons/1e-8/DMaxO;
     
     sprintf( F,"%s/.D_%s_%02d",DIRECTORY,HSTFile,NProj );
     OF.open( F );
     if(!OF)
      {cerr << "# Error: SaveHistoR, Can't create "<<F<<"\n";exit(WRITE_ERROR);}
      
     OF << "#  double ions = NumberOfImplants*EffectiveWeigth, Virtual ions\n";
     OF << "# " << RealIons << " " << LVI.NumberOfElements() << endl;
     OF << "#  Distance_travelled(nm), Conc. (at/cm^3\n";
     for(i=0; i<N_BIN; i++ )
      OF << i*DMax/N_BIN*0.1 << " " << HR[i]*K << endl;
     OF.close(); 
   }

void Hist::ShowFinalInfo(int NProj)
  {
     std::cout << endl;
     std::cout << "Projectile " << NProj
          << ",\tdouble ions = " << RealIons 
          << ",\tVirtual ions = " << LVI.NumberOfElements() << endl;
  }

void Hist::SolvePearsonIV(int NProj,char * HSTFile)
  {      
   // File format: 
   //   depth (nm), conc.Gauss (at/cm^2), conc.PearsonIV (at/cm^2) 
   ofstream ES;
   int i;
   double m1,m2,m3,m4, G[N_BIN+1], P[N_BIN+1], xp;
   double x, inte=0, inte2 = 0;

   // Solve the four moments of the Pearson IV distribution
   m1=m2=m3=m4=0.0;     
   for( i=0; i< N_BIN; i++ ) m1 = m1 + H[i]*i;
   m1 = m1 * DMax/N_BIN/RealIons;

   for( i=0; i<N_BIN; i++ )
   {
    xp  = i*DMax/N_BIN - m1;
    m2  = m2 + xp*xp * H[i];
    m3  = m3 + xp*xp*xp * H[i];
    m4  = m4 + xp*xp*xp*xp * H[i]; 
   }    
   m2 = m2 / RealIons;
   m3 = m3 / RealIons;
   m4 = m4 / RealIons;
   m2=sqrt(m2);
   if( m2 != 0)
    { 
     m3 = m3 / (m2*m2*m2);
     m4 = m4 / (m2*m2*m2*m2);
    }
   char F[300];
   sprintf(F,"%s/PearsonIV_%02d.dat",DIRECTORY,NProj);
   ES.open( F );
   ES << "#  Pearson IV moments:" <<  endl;
   ES << "# Mean range m1 = " << m1 << endl;
   ES << "# Straggle   m2 = " << m2 << endl;
   ES << "# Skewness   m3 = " << m3 << endl;
   ES << "# Kurtosis   m4 = " << m4 << endl;
  
   double A  = 10*m4-12*m3*m3-18;
   double a  = - m3*m2*(m4+3)/A;
   double b0 = - m2*m2*(4*m4-3*m3*m3)/A;
   double b1 = a;
   double b2 = -(2*m4-3*m3*m3-6)/A;

   if(  b1*b1-4*b0*b2 >=0 )
   {
    ES << "# Cannot generate a Pearson IV" << endl;
    return;
   }
   for( x = 0,i=0; x < DMax; x+= DMax/N_BIN,i++ )
    { 
     xp = x-m1 ;
     double s1 = log( fabs(b2*xp*xp+b1*xp+b0)/m1/m1 )/(2*b2);
     double s2 = (b1/b2+2*a)/sqrt(4*b2*b0-b1*b1);
     double s3 = atan2(2*b2*xp+b1,sqrt(4*b2*b0-b1*b1));
     double nn = exp( (s1 - s2*s3) );         
     inte += nn*DMax/N_BIN*1e-8;
     double gaus = 1/sqrt(2.0*M_PI)/m2*exp( -(x-m1)*(x-m1)/2/m2/m2 );
     inte2+= gaus*DMax/N_BIN*1e-8;
     G[i]=gaus;
     P[i]=nn;
    }
   for( x = 0, i=0; x < DMax; x+= DMax/N_BIN,i++ )
    { 
     ES << x*0.1 << " " 
        << Dose * G[i] / inte2 << " "
        << Dose * P[i] / inte  << endl;
    }
   ES.close();
  }

//=========================================================================

// Add a virtual ion to the list
void Histo3D::Add( Vector R, double W, double O )
 {
  IonVirtual IV;

  IV.R      = R;
  IV.Weigth = W;   
  IV.Other  = O;
  LIV.Add(IV);
 }

void Histo3D::Show( int NProj, double Dose )
{
  LinkObject<IonVirtual> * IV;
  FILE *OF;
  double RealIons = 0.0;
      
  char F[300];
  sprintf( F,"%s/Ion.%02d.coordinates.dat",DIRECTORY, NProj );
  OF= fopen( F,"w" );
  if(OF)
   {
    fprintf(OF,"# Format\n");
    fprintf(OF,"#    X (A)      Y (A)      Z (A) Stat. Weigth Dist. travel. (A)\n");
    IV=LIV.PointerToFirstData();
    do
     {
      fprintf(OF,"%10.2lf %10.2lf %10.2lf %12.6lf %10.2lf\n",
            IV->Data.R.X, IV->Data.R.Y, IV->Data.R.Z, IV->Data.Weigth, IV->Data.Other);
      RealIons += IV->Data.Weigth;
      IV=LIV.PointerToNextData();
     } while(IV != NULL);
    fprintf(OF,"# Dose (at/cm^2) Real_Ions Virtual_Ions\n");
    fprintf(OF,"# %14g %9.0lf %12d\n", Dose, RealIons, LIV.NumberOfElements());
    fclose(OF);
   }
}

// Process the virtual ion list to solve the histogram (3D)
void Histo3D::Solve3D( int NProj, double Dose, double WinAB, int ascii )
 {
  // << Initialize >>
  int i,j,k;
  for(i=0;i<N_BIN2D;i++)
   for(j=0;j<N_BIN2D;j++)
    for(k=0;k<N_BIN2D;k++) H3D[i][j][k] = 0.0;
   
  // << Find the extremes >>
  LinkObject<IonVirtual> * IV;
  double XMax = 0.0;
  double YMax = 0.0;
  double YMin = 0.0;
  double ZMax = 0.0;
  double ZMin = 0.0;
  IV=LIV.PointerToFirstData();
  do
   {
     Vector R = OR.Proyecta( IV->Data.R );

     if( R.X>XMax ) XMax = R.X;
     if( R.Y>YMax ) YMax = R.Y;
     if( R.Y<YMin ) YMin = R.Y;
     if( R.Z>ZMax ) ZMax = R.Z;
     if( R.Z<ZMin ) ZMin = R.Z;

     IV=LIV.PointerToNextData();
   }
  while(IV != NULL);

  // << Solve >>
  double Kx =  XMax      /N_BIN2D; 
  double Ky = (YMax-YMin)/N_BIN2D;
  double Kz = (ZMax-ZMin)/N_BIN2D;
  double RealIons = 0.0;
  double OutIons  = 0.0;
  IV=LIV.PointerToFirstData();
  do
   {
     Vector R = OR.Proyecta( IV->Data.R );

     i = (int) ( R.X       / Kx );
     j = (int) ((R.Y-YMin) / Ky );
     k = (int) ((R.Z-ZMin) / Kz );
     
     if(R.X>=0) H3D[i][j][k] +=          IV->Data.Weigth;
     else        OutIons      = OutIons + IV->Data.Weigth;
     RealIons  += IV->Data.Weigth;

     IV=LIV.PointerToNextData();
   }
  while(IV != NULL);

  // << Save the histogram >> 
  if( RealIons <=0 ) { cerr << "RealIons <=0\n"; exit(HISTOGRAM_ERROR); };

  double K  = Dose/RealIons * WinAB/(Kx*Ky*Kz)*1e8;

  FILE *OF, *OF2;
  char F[300];
  if(ascii)     // Optional ASCII output =================
  {
   sprintf( F,"%s/Ion.%02d.histo3D.dat",DIRECTORY, NProj );
   OF = fopen( F,"w" );
   if(OF!=NULL)
   {
    fprintf(OF,"# double Ions    = %.1lf\n", RealIons);
    fprintf(OF,"# Virtual Ions = %d\n",  LIV.NumberOfElements());
    OutIons =OutIons*100/RealIons;
    fprintf(OF,"# double Ions outside = %4.2lf %%\n", OutIons );
    fprintf(OF,"# Volume cell  = %lg cm^3\n", Kx*Ky*Kz*1e-24);
    fprintf(OF,"# Format:\n");
    fprintf(OF,"#   X (nm)     Y (nm)     Z (nm) Conc(at/cm^3)\n");
    for( i=0; i<N_BIN2D; i++ )
     {
      for( j=0; j<N_BIN2D; j++ )
       { 
        for( k=0; k<N_BIN2D; k++ )
         {
          fprintf(OF,"%10.2lf %10.2lf %10.2lf %13.4lg\n",
           Kx*i*0.1, (Ky*j+YMin)*0.1, (Kz*k+ZMin)*0.1, K*H3D[i][j][k]);
         }
        fprintf(OF,"\n"); 
       }
      fprintf(OF,"\n"); 
     }
    fclose(OF);
   }
  }
  // Binary output =======================================
  sprintf( F,"%s/Ion.%02d.histo3D.bin",DIRECTORY, NProj );
  OF2 = fopen( F,"wb" );
  if(OF2!=NULL)
   {
    for( i=0; i<N_BIN2D; i++ )
      for( j=0; j<N_BIN2D; j++ )
        for( k=0; k<N_BIN2D; k++ )
         {          
          H3D[i][j][k] = K*H3D[i][j][k];
         }         
    int auxi;
    double auxf;     
    auxi = N_BIN2D; fwrite( &auxi, sizeof(int   ), 1, OF2); // N_X
    auxi = N_BIN2D; fwrite( &auxi, sizeof(int   ), 1, OF2); // N_Y
    auxi = N_BIN2D; fwrite( &auxi, sizeof(int   ), 1, OF2); // N_Z
    auxf = 0.0;     fwrite( &auxf, sizeof(double), 1, OF2); // X
    auxf = XMax;    fwrite( &auxf, sizeof(double), 1, OF2);
    auxf = YMin;    fwrite( &auxf, sizeof(double), 1, OF2); // Y
    auxf = YMax;    fwrite( &auxf, sizeof(double), 1, OF2);
    auxf = ZMin;    fwrite( &auxf, sizeof(double), 1, OF2); // Z
    auxf = ZMax;    fwrite( &auxf, sizeof(double), 1, OF2);

    fwrite(H3D,sizeof(double ),N_BIN2D*N_BIN2D*N_BIN2D,OF2);
    fclose(OF2);
   }     
 }

 Hist     HstPROJ[N_ATOMS];    // Statistics 1D
 Histo3D  HstPROJ3D[N_ATOMS];  // Statistics 3D

}



