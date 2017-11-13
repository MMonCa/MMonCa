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
////////////////////////////////////////////////////////////////////////
//
//  AMORF3D.CPP
//  Modified Kinchin-Pease model for Damage Accumulation
//
////////////////////////////////////////////////////////////////////////
#include "amorf3d.h"

using std::cerr;
using std::cout;
using std::ifstream;
using std::endl;

namespace BCA {
//  Constructor ---------------------------------------------
clase_Amorfizacion3D::clase_Amorfizacion3D()
  {
   DefectDensity.Set3D();
   DefectDensity.Clear3D();
#ifdef FAST_KP
   E_Acc.Set3D(); E_Acc.Clear3D(0.0);
#else   
   AA.Set3D();   AA.Clear3D(1.0);
   BB.Set3D();   BB.Clear3D(0.0);
#endif
   k    =   0.8;        // Kinchin-Pease constant
   frec =   0.06;       // Surviving recombination factor
   Na   =   0.0049968;  // Amorphization density (for Si)
   NNN  =   0;          // Number of damage accumulation operations

   min_y = min_z = N_SECTOR;
   max_y = max_z = 0;
  
#ifdef ENERGIA_DEPOSITADA
   ETOTAL = 0.0;
   OO.open("Energia_Depositada");
#endif
#ifdef MODEL2
   DeltaDamageByCascade=0;
#endif
 }

// Init -------------------------------------------------------------------
void clase_Amorfizacion3D::Init3D( double Dose, double IA, int LS,
        double RecombinationFactor, double KinchinPeaseConstant,
        double Na_, double COE, double DE, double WA, double WB,
        int MultipleImplant, const std::string &PreviousDamageFile)
  {
   ifstream In;
   double depthx,depthy,depthz, percentage, lado;

   frec              = RecombinationFactor;
   k                 = KinchinPeaseConstant;
   Na                = Na_;
   CutOffEnergy      = COE;
   LevelOfSimulation = LS;
   DisplacementEnergy= DE;

   Area = IA;  // (A^2)  
   lado = sqrt(Area);

   if(lado > WA) WinA = lado; else WinA = WA;
   if(lado > WB) WinB = lado; else WinB = WB;
    
   DefectDensity.mExtremoY = lado*(N_SECTOR/2); DefectDensity.mAnchoY = lado;
   DefectDensity.mExtremoZ = lado*(N_SECTOR/2); DefectDensity.mAnchoZ = lado;   
#ifdef FAST_KP
   KK3  = k/(2*DisplacementEnergy) * frec;
   E_Acc.mExtremoY = lado*(N_SECTOR/2);   E_Acc.mAnchoY   = lado;
   E_Acc.mExtremoZ = lado*(N_SECTOR/2);   E_Acc.mAnchoZ   = lado;
#else
   KK2  = k/(2*DisplacementEnergy) * frec/ Area;
   KK   = KK2/Na;
   AA.mExtremoY = lado*(N_SECTOR/2);   AA.mAnchoY   = lado;
   AA.mExtremoZ = lado*(N_SECTOR/2);   AA.mAnchoZ   = lado;
   BB.mExtremoY = lado*(N_SECTOR/2);   BB.mAnchoY   = lado;
   BB.mExtremoZ = lado*(N_SECTOR/2);   BB.mAnchoZ   = lado;
#endif
   // Multiple implant inicialization -------------------------------
   if( MultipleImplant==0 )
   {
#ifdef FAST_KP
    E_Acc.Clear3D(0.0);    
#else   
    AA.Clear3D(1.0);
    BB.Clear3D(0.0);
#endif    
   }
   // Read previous damage if exits ---------------------------------
   if (PreviousDamageFile.size())
    {
      In.open(PreviousDamageFile.c_str());
      if(In == NULL)
        {
         cerr << "# I can't open ("<<PreviousDamageFile<<")\n";
         exit(READ_ERROR);
        }
      cout << "# Reading ("<<PreviousDamageFile<<")\n";
      int i,j,k;
      while(!In.eof())
      for(i=0;i<N_BINS;i++)
       for(j=0;j<N_SECTOR;j++)
        for(k=0;k<N_SECTOR;k++)
        {
         In >> depthx >> depthy >> depthz 
            >> DefectDensity.M[i][j][k] >> percentage;
        }
      DefectDensity.Set3D(0.0,depthx,depthy,depthz);
    }
  }

 // Damage accumulation ---------------------------------------------------
#ifdef FAST_KP
void clase_Amorfizacion3D::DamageAccumulation3D
       ( Vector pos, double energia, char tipo, double Weigth )
  {
   Vector l;
   l = E_Acc.GetIndex( pos );
   double Dens_energia = energia/DefectDensity.mAncho;
   if(l.X>=0)
    {
     if(LevelOfSimulation==1)
      { // Simple case. Follow only the primary ion =============
        if(Dens_energia>CutOffEnergy*Area) 
                        energia = CutOffEnergy*Area*DefectDensity.mAncho;
      }
     else
      { // Complex case. Complete cascade =======================
        if(tipo=='B') energia=DisplacementEnergy*Weigth;
      }
     E_Acc.M[(int)l.X][(int)l.Y][(int)l.Z] = 
     E_Acc.M[(int)l.X][(int)l.Y][(int)l.Z] + energia;

     if( l.Y > max_y) max_y = l.Y;
     if( l.Z > max_z) max_z = l.Z;
     if( l.Y < min_y) min_y = l.Y;
     if( l.Z < min_z) min_z = l.Z;

    }
  };
  
#else  
void clase_Amorfizacion3D::DamageAccumulation3D
       ( Vector pos, double energia, char tipo ,double Weigth )
  {
   Vector l;
   l = AA.GetIndex( pos );
   l = BB.GetIndex( pos );
   if(l.X>=0)
    {
     energia=energia/AA.mAncho;
     if(LevelOfSimulation==1)
      { // Simple case. Follow only the primary ion =============
        if(energia>CutOffEnergy*Area) energia = CutOffEnergy*Area;
      }
     else
      { // Complex case. Complete cascade =======================
        if(tipo=='B') energia=DisplacementEnergy*Weigth/AA.mAncho;
      }
     AA.M[(int)l.X][(int)l.Y][(int)l.Z]=
     AA.M[(int)l.X][(int)l.Y][(int)l.Z]*(1-KK*energia);
     BB.M[(int)l.X][(int)l.Y][(int)l.Z]=
     BB.M[(int)l.X][(int)l.Y][(int)l.Z]*(1-KK*energia)+ KK2*energia;

     if( l.Y > max_y) max_y = l.Y;
     if( l.Z > max_z) max_z = l.Z;
     if( l.Y < min_y) min_y = l.Y;
     if( l.Z < min_z) min_z = l.Z;
    }
  }
#endif

 // Update defect density -------------------------------------------------
#ifdef FAST_KP
void clase_Amorfizacion3D::UpdateDefectDensity3D(Vector InitPos)
  {
   int   jj,kk;
   double y,z;
   int i,j,k;
   
   DefectDensity.ReScaleX(0.0,E_Acc.mUltimo);
   double K = KK3/Area/DefectDensity.mAncho;

   double min_InitPosY = -WinA/2-InitPos.Y;
   double max_InitPosY =  WinA/2-InitPos.Y;
   double min_InitPosZ = -WinB/2-InitPos.Z;
   double max_InitPosZ =  WinB/2-InitPos.Z;
      
//1�.- Scan all zones
   y = -E_Acc.mExtremoY;
     
   for(jj=-N_SECTOR/2;jj<N_SECTOR/2;jj++)
   {
    z = -E_Acc.mExtremoZ;

    for(kk=-N_SECTOR/2;kk<N_SECTOR/2;kk++)
    {
      //2�.- Zone validation (**) Care with bit zones
      if(y>=min_InitPosY && y<=max_InitPosY && 
         z>=min_InitPosZ && z<=max_InitPosZ)    
      {       
        //3�.- Update defect density in this zone
        for(j=0;j<N_SECTOR;j++)
         for(k=0;k<N_SECTOR;k++)
          {
            int la_j= j+jj;
            int la_k= k+kk;

            if(la_j>=min_y && la_j<=max_y && la_k>=min_z && la_k<=max_z )           
             for(i=0;i<N_BINS;i++)
             {
              DefectDensity.M[i][j][k] +=  
                  K*E_Acc.M[i][la_j][la_k]*( 1 - DefectDensity.M[i][j][k]/Na );
             }
           }      
      }
      z += E_Acc.mAnchoZ;
    }
    y += E_Acc.mAnchoY; 
   } 
//4�.- Matrix initialization
/*
   for(i=0;i<N_BINS;i++)
    for(j=0;j<N_SECTOR;j++)
     for(k=0;k<N_SECTOR;k++) E_Acc.M[i][j][k]=0.0;
*/
   t_double * PP;
   PP = &E_Acc.M[0][0][0];
   for(i=0;i<N_BINS*N_SECTOR*N_SECTOR;i++) // Faster
     *PP++=0.0;

   if(NNN%100==0 && NNN!=0) WriteAmorphization3D(0);
   NNN++;
  }
#else
void clase_Amorfizacion3D::UpdateDefectDensity3D(Vector InitPos)
  {
   double y,z;
   int jj,kk;   
   int i,j,k;
   DefectDensity.ReScaleX(0.0,AA.mUltimo);
 #ifndef DAMAGE3D_FASE1

   double min_InitPosY = -WinA/2-InitPos.Y;
   double max_InitPosY =  WinA/2-InitPos.Y;
   double min_InitPosZ = -WinB/2-InitPos.Z;
   double max_InitPosZ =  WinB/2-InitPos.Z;
      
//1�.- Scaning all th zones
   y = -AA.mExtremoY;
   for(jj=-N_SECTOR/2;jj<N_SECTOR/2;jj++)
   {
    z = -AA.mExtremoZ;
    for(kk=-N_SECTOR/2;kk<N_SECTOR/2;kk++)
    {
      //2�.- Validate the zone (**) Care with bit zones
      if(y>=min_InitPosY && y<=max_InitPosY && 
         z>=min_InitPosZ && z<=max_InitPosZ)    
      {       
        //3�.- Update defect density for that zone
        for(j=0;j<N_SECTOR;j++)
         for(k=0;k<N_SECTOR;k++)
          {
            int la_j= j+jj;
            int la_k= k+kk;
            if(la_j>=min_y && la_j<=max_y && la_k>=min_z && la_k<=max_z )            
             for(i=0;i<N_BINS;i++)
             {
              DefectDensity.M[i][j][k] = 
                AA.M[i][la_j][la_k]*DefectDensity.M[i][j][k] 
                + BB.M[i][la_j][la_k];
             }
           }      
      }
      z += AA.mAnchoZ;
    }
    y += AA.mAnchoY; 
   } 
//4�.- Matrix initialization
   for(i=0;i<N_BINS;i++)
    for(j=0;j<N_SECTOR;j++)
     for(k=0;k<N_SECTOR;k++)
     {
      AA.M[i][j][k] = 1.0;
      BB.M[i][j][k] = 0.0;
     }
 #else
   for(i=0;i<N_BINS;i++)
    for(j=0;j<N_SECTOR;j++)
     for(k=0;k<N_SECTOR;k++)
    {
     DefectDensity.M[i][j][k] = AA.M[i][j][k]*DefectDensity.M[i][j][k] + BB.M[i][j][k];
     AA.M[i][j][k] = 1.0;
     BB.M[i][j][k] = 0.0;
    }
 #endif  
   if(NNN%100==0 && NNN!=0) WriteAmorphization3D(0);
   NNN++;
  }
#endif

 // Damage accumulation. Sequential version -------------------------------
void clase_Amorfizacion3D::SeqDamageAccumulation3D
       ( Vector pos, double energia, char tipo ,double Weigth, double Temp )
  {
   double delta_n, n, volume;

   Vector l=DefectDensity.GetIndex(pos);
   if(l.X>=0)
    {
     volume = Area*DefectDensity.mAncho;
     if(LevelOfSimulation==1)
      { if( energia > CutOffEnergy*volume ) energia = CutOffEnergy*volume; }
     else
      { if( tipo=='B')               energia = DisplacementEnergy*Weigth; }
     n       = k * energia/( 2*DisplacementEnergy );
     
     #ifdef MODEL2
     delta_n = n * frec * (1-DefectDensity.M[(int)l.X][(int)l.Y][(int)l.Z]/Na);     
     DeltaDamageByCascade+=delta_n;
     #else
     delta_n = n * frec * (1-DefectDensity.M[(int)l.X][(int)l.Y][(int)l.Z]/Na);
     #endif
     
     DefectDensity.M[(int)l.X][(int)l.Y][(int)l.Z]
         = DefectDensity.M[(int)l.X][(int)l.Y][(int)l.Z]+delta_n/volume;
    }
   if(NNN%10000==0) WriteAmorphization3D(0);
   NNN++;
 }

 // Write amorphization histogram -----------------------------------------
void clase_Amorfizacion3D::WriteAmorphization3D(int Opciones)
{
   FILE * OF;
   int i,j,k;
   char F[300];

   // Output 3D Damage in Binary ================================

   sprintf(F, "%s/Target.damage3D.bin", DIRECTORY );
   OF= fopen(F,"wb");
   if(!OF) { cerr << "# I can't write in " << F << endl; exit(WRITE_ERROR); };

   int auxi;
   double auxf;     
   auxi = N_BINS;   fwrite( &auxi, sizeof(int   ), 1, OF); // N_X
   auxi = N_SECTOR; fwrite( &auxi, sizeof(int   ), 1, OF); // N_Y
   auxi = N_SECTOR; fwrite( &auxi, sizeof(int   ), 1, OF); // N_Z
   auxf = 0.0;      fwrite( &auxf, sizeof(double), 1, OF); // X
   auxf =      N_BINS*DefectDensity.mAncho;  fwrite( &auxf, sizeof(double), 1, OF);
   auxf = -N_SECTOR/2*DefectDensity.mAnchoY; fwrite( &auxf, sizeof(double), 1, OF); // Y
   auxf =  N_SECTOR/2*DefectDensity.mAnchoY; fwrite( &auxf, sizeof(double), 1, OF);
   auxf = -N_SECTOR/2*DefectDensity.mAnchoZ; fwrite( &auxf, sizeof(double), 1, OF); // Z
   auxf =  N_SECTOR/2*DefectDensity.mAnchoZ; fwrite( &auxf, sizeof(double), 1, OF);

   for(i=0;i<N_BINS;i++)
    for(j=0;j<N_SECTOR;j++)
     for(k=0;k<N_SECTOR;k++)
     {
      double realTemp = DefectDensity.M[i][j][k]*100/Na;
      fwrite( &realTemp, sizeof(double),1, OF );
     }
//   fwrite( DefectDensity.M, sizeof(double),N_BINS*N_SECTOR*N_SECTOR, OF );
   fclose(OF);

   // Output 3D damage file in ASCII ================================
 if(Opciones)
 {
   sprintf(F, "%s/Target.damage3D.dat", DIRECTORY );
   OF= fopen(F,"w");
   if(!OF) { cerr << "# I can't write in " << F << endl; exit(WRITE_ERROR); };
   fprintf(OF,"# 3D damage file %d x %d x %d boxes\n",N_BINS,N_SECTOR,N_SECTOR);
   fprintf(OF,"#   X (nm)      Y(nm)      Z(nm) Defects_density (A^-3) %%Amorphization\n");
   for(i=0;i<N_BINS;i++)
     for(j=0;j<N_SECTOR;j++)
      {
       for(k=0;k<N_SECTOR;k++)
        fprintf(OF,"%10.2lf %10.2lf %10.2lf %12.2lf %18.2lf\n", 
                         i*DefectDensity.mAncho /10,
            (j-N_SECTOR/2)*DefectDensity.mAnchoY/10,
            (k-N_SECTOR/2)*DefectDensity.mAnchoZ/10,
                           DefectDensity.M[i][j][k], 
                    DefectDensity.M[i][j][k]*100/Na);
       fprintf(OF,"\n");
      };
   fclose(OF);
 }
}

 // SUMA Energia ----------------------------------------------------------
#ifdef ENERGIA_DEPOSITADA
void clase_Amorfizacion3D::Suma(double X, double Y, double Z, double Energy,
               char C[2], double Weigth)
  {
   if(strcmp(C,"B")==0)  ETOTAL = ETOTAL + DisplacementEnergy*Weigth;
   else                  ETOTAL = ETOTAL + Energy;
   OO << X << " " << Y << " " << Z << " " << Energy << " " << C << endl;
  };
#endif  

 // Destructor ------------------------------------------------------------
clase_Amorfizacion3D::~clase_Amorfizacion3D()
   {
#ifdef ENERGIA_DEPOSITADA
    cout << "####### Total Energy = " << ETOTAL << endl;
#endif
//    WriteAmorphization3D(0);
   }

clase_Amorfizacion3D Amorfizacion3D;
}
// ------------------------------------------------------------------------
