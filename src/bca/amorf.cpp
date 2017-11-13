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
//  AMORF.CPP
//  Modified Kinchin-Pease model for Damage Accumulation
//
////////////////////////////////////////////////////////////////////////
#include "amorf.h"

using std::ifstream;
using std::ofstream;
using std::cerr;
using std::endl;
using std::cout;

namespace BCA {

//  Constructor ---------------------------------------------
clase_Amorfizacion::clase_Amorfizacion()
  {
   DefectDensity.Set();
   DefectDensity.Clear();
#ifdef FAST_KP
   E_Acc.Set(); E_Acc.Clear(0.0);
#else   
   AA.Set();   AA.Clear(1.0);
   BB.Set();   BB.Clear(0.0);
#endif
   k    =   0.8;        // Kinchin-Pease constant
   frec =   0.06;       // Surviving recombination factor
   Na   =   0.0049968;  // Amorphization density (for Si)
   NNN  =   0;          // Number of damage accumulation operations

#ifdef MODEL2
   DeltaDamageByCascade = 0;
   CteB = 0.43;
   CteD = 0.1;
   CteK = 0.1;
   Mean_CrossSection = 0; // cm^2
   
#endif
#ifdef ENERGIA_DEPOSITADA
   ETOTAL = 0.0;
   OO.open("Distributed_energy.dat");
#endif 
 }

// Init -------------------------------------------------------------------
void clase_Amorfizacion::Init( double Dose, double IA, int LS,
        double RecombinationFactor, double KinchinPeaseConstant,
        double Na_, double COE, double DE,
        int MultipleImplant, const std::string &PreviousDamageFile, double DR, double WE )
  {
   ifstream In;
   double depth, percentage;

#ifdef MODEL2
   if(ancho>0.0) DefectDensity.Set(0.0, ancho*N_BINS);
   Weff = WE;
#endif

   frec              = RecombinationFactor;
   k                 = KinchinPeaseConstant;
   Na                = Na_;
   CutOffEnergy      = COE;
   LevelOfSimulation = LS;
   DisplacementEnergy= DE;

   DoseRate          = DR;

   Area = IA;                                   // (A^2)
#ifdef FAST_KP
   KK3  = k/(2*DisplacementEnergy) * frec;
#else
   KK2  = k/(2*DisplacementEnergy) * frec/ Area;
   KK   = KK2/Na;
#endif
   // Multiple implant inicialization -------------------------------
   if( MultipleImplant==0 )
   {
#ifdef FAST_KP
    E_Acc.Clear(0.0);    
#else   
    AA.Clear(1.0);
    BB.Clear(0.0);
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
      int i=0; char C;
      In >> C 
         >> DefectDensity.mPrimero 
         >> DefectDensity.mUltimo 
         >> DefectDensity.mAncho;
      for(i=0; i<N_BINS; i++)
         In >> depth >> DefectDensity.M[i] >> percentage;
      In.close();
      #ifdef FAST_KP
      E_Acc.mPrimero=DefectDensity.mPrimero;    
      E_Acc.mUltimo =DefectDensity.mUltimo;
      E_Acc.mAncho  =DefectDensity.mAncho;
      #else   
      AA.mPrimero=DefectDensity.mPrimero;    
      AA.mUltimo =DefectDensity.mUltimo;
      AA.mAncho  =DefectDensity.mAncho;
      BB.mPrimero=DefectDensity.mPrimero;    
      BB.mUltimo =DefectDensity.mUltimo;
      BB.mAncho  =DefectDensity.mAncho;
      #endif    
   }
  }

 // Damage accumulation ---------------------------------------------------
#ifdef FAST_KP
void clase_Amorfizacion::DamageAccumulation
       ( double posx, double energia, char tipo, double Weigth )
  {
   int l;
   l = E_Acc.GetIndex( posx );
   double Dens_energia = energia/E_Acc.mAncho;
   if(l>=0)
    {
     if(LevelOfSimulation==1)
      { // Simple case. Follow only the primary ion =============
        if(Dens_energia>CutOffEnergy*Area) 
                        energia = CutOffEnergy*Area*E_Acc.mAncho;
      }
     else
      { // Complex case. Complete cascade =======================
        if(tipo=='B') energia=DisplacementEnergy*Weigth;
      }
     E_Acc.M[l] = E_Acc.M[l] + energia;
    }
  };
  
#else  
void clase_Amorfizacion::DamageAccumulation
       ( double posx, double energia, char tipo ,double Weigth )
  {
   int l;
   l = AA.GetIndex( posx );
   l = BB.GetIndex( posx );
   if(l>=0)
    {
     energia=energia/AA.mAncho;
     if(LevelOfSimulation==1)
      { // Simple case. Follow only the primary ion =============
        if(energia>CutOffEnergy*Area) 
                        energia = CutOffEnergy*Area;
      }
     else
      { // Complex case. Complete cascade =======================
        if(tipo=='B') energia=DisplacementEnergy*Weigth/AA.mAncho;
      }
     AA.M[l]=AA.M[l]*(1-KK*energia);
     BB.M[l]=BB.M[l]*(1-KK*energia)+ KK2*energia;     
    }
  }
#endif

 // Update defect density -------------------------------------------------
#ifdef FAST_KP
void clase_Amorfizacion::UpdateDefectDensity()
  {
   int l;
   DefectDensity.ReScale(0.0,E_Acc.mUltimo);
   double K = KK3/Area/DefectDensity.mAncho;
   for(l=0;l<N_BINS;l++)
    {
     DefectDensity.M[l] = DefectDensity.M[l] 
                            +  K*E_Acc.M[l]*( 1 - DefectDensity.M[l]/Na );
     E_Acc.M[l]=0.0;
    }
   if(NNN%100==0 && NNN!=0) WriteAmorphization();
   NNN++;
  }
#else
void clase_Amorfizacion::UpdateDefectDensity()
  {
   int l;
   DefectDensity.ReScale(0.0,AA.mUltimo);
   for(l=0;l<N_BINS;l++)
    {
    #ifdef LOURDES_11_02
     double Old_DD= DefectDensity.M[l];
    #endif
     DefectDensity.M[l] = AA.M[l]*DefectDensity.M[l] + BB.M[l];
     AA.M[l] = 1.0;
     BB.M[l] = 0.0;
    #ifdef LOURDES_11_02
     DeltaN += (DefectDensity.M[l]-Old_DD)*Area*DefectDensity.mAncho;
    #endif
    }
   if(NNN%100==0 && NNN!=0) WriteAmorphization();
   NNN++;
  }
#endif

 // Damage accumulation. Sequential version -------------------------------
void clase_Amorfizacion::SeqDamageAccumulation
       ( double posx, double energia, char tipo ,double Weigth, double Temp )
  {
   double delta_n, n, volume;

   int l=DefectDensity.GetIndex(posx);
   if(l>=0)
    {
     volume = Area*DefectDensity.mAncho;
     if(LevelOfSimulation==1)
      { if( energia > CutOffEnergy*volume ) energia = CutOffEnergy*volume; }
     else
      { if( tipo=='B')               energia = DisplacementEnergy*Weigth; }
     n       = k * energia/( 2*DisplacementEnergy );
    
    
     #ifdef MODEL2
/*
FILE *OO;
OO=fopen("E_mean","a+");
fprintf(OO,"%d %g\n",l,energia);
fclose(OO);
*/

     
     TotalGeneratedFrenkelPairs += n;
     
     double kB = 1.38e-23; // Julios/K
     
//     frec = CteA*(1.0 -exp(-CteB / (kB*Temp) ));

// Model M
//     frec = 1/(1.0 +exp(CteA*(Temp-CteB) ));

// Model N
//     frec =CteA/(1.0 +exp(Temp-CteB));

// Modelo O
//     frec = CteA*( 1-exp(-(DefectDensity.M[l]+n/volume)/Na) )/(1.0 +exp(Temp-CteB));     
//     delta_n = frec * n * (1-DefectDensity.M[l]/Na);
     
// Modelo Q
//     frec = ( 1-exp(-CteA*(DefectDensity.M[l]+n/volume)/Na) )/(1.0 +exp( 0.1*(Temp-CteB))); 
//     delta_n = frec * n * (1-DefectDensity.M[l]/Na);
     
    
// Modelo R    
//       double Tc = 1/(CteB+CteC*log10(DoseRate));
//       frec = ( 1-exp(-CteA*(DefectDensity.M[l]+n/volume)/Na) )/(1.0 +exp( 0.1*(Temp-Tc))); 
//       delta_n = frec * n * (1-DefectDensity.M[l]/Na);

// Modelo R'
//       double Tc = 1/(CteB+CteC*log10(DoseRate));
//       n = n / Weff;       
//       frec = ( 1-exp(-CteA*(DefectDensity.M[l]+n/volume)/Na) )/(1.0 +exp( 0.1*(Temp-Tc))); 
//       delta_n = frec * n * (1-DefectDensity.M[l]/Na);

     

// Modelo R2   
//       double Tc = 1/(CteB+CteC*log10(DoseRate));
//       frec = ( 1-exp(-CteA*(DefectDensity.M[l]+n/volume)/Na) )/(1.0 +exp( 0.1*(Temp-Tc))); 
//       delta_n = frec * n * (1-DefectDensity.M[l]/Na);

// Modelo R3
//       double Tc = 1/(CteB+CteC*log10(DoseRate));
//       frec = CteA/(1.0 +exp( 0.1*(Temp-Tc))); 
//       delta_n = frec * n * (1-DefectDensity.M[l]/Na);

// Modelo R4
//       double Tc = 1/(CteB+CteC*log10(DoseRate));
//       frec = ( 0.2+0.15*exp( 4.0 *DefectDensity.M[l]/Na) )/(1.0 +exp( 0.1*(Temp-Tc))); 
//       delta_n = frec * n * (1-DefectDensity.M[l]/Na);

// Modelo R5

//       double Tc = 1/(CteB+CteC*log10(DoseRate));
//       double E = (1-CteD)/(1-exp(-CteA));
//       frec    =( CteD + E*( 1-exp(-CteA*DefectDensity.M[l]/Na) ))
//    		/(1.0 +exp( CteK*(Temp-Tc))); 
//       delta_n = frec * n * (1-DefectDensity.M[l]/Na);


// Modelo T    
//     double Tc = 1/(CteB+CteC*log10(DoseRate));
//     double A = 1.60+2.628*log(DoseRate/1.0e11);
//     frec = ( 1 - exp( - A*(DefectDensity.M[l]+n/volume)/Na) )/(1.0 +exp( 0.1*(Temp-Tc))); 
//     delta_n = frec * n * (1-DefectDensity.M[l]/Na);
     

// Modelo U
//     double Tc = 1/(CteB+CteC*log10(DoseRate));
//     frec = ( 1 - exp( - CteA/(kB*Temp)*(DefectDensity.M[l]+n/volume)/Na) )/(1.0 +exp( 0.1*(Temp-Tc))); 
//     delta_n = frec * n * (1-DefectDensity.M[l]/Na);
     

// Modelo V
//     double Tc = 1/(CteB+CteC*log10(DoseRate));
//     double A = 20-0.0825*(Temp-300.0);
//     frec = ( 1 - exp( - A*(DefectDensity.M[l]+n/volume)/Na) )/(1.0 +exp( 0.1*(Temp-Tc))); 
//     delta_n = frec * n * (1-DefectDensity.M[l]/Na);
     

// Modelo W
/*  
       //                t_0  * exp( E_A              / (k  * T   ))
       double t_rec = 162e-15 * exp( 0.43 * 1.602e-19 / (kB * Temp));

//       frec    = CteA * 
//             DoseRate * Mean_CrossSection * t_rec;
       frec    = CteD *  DoseRate  * t_rec;
       delta_n = frec * n * (1-DefectDensity.M[l]/Na);
*/
// Modelo W2
//     double rN = DefectDensity.M[l]/Na;

//     double t_rec = 162e-15 * exp( (CteB+CteC*rN) * 1.602e-19 / (kB * Temp));
//     frec    = CteD *  DoseRate  * t_rec;

//     if(frec>1.0) frec=1.0;
//     delta_n = frec * n * (1-rN);

// Modelo W3
     double rN = DefectDensity.M[l]/Na;

     double t_rec = (162e-15+CteC*rN) * exp( (CteB) * 1.602e-19 / (kB * Temp));
     frec    = CteD *  DoseRate  * t_rec;

     if(frec>1.0) frec=1.0;
     delta_n = frec * n * (1-rN);


       
     DeltaDamageByCascade += delta_n;
     #else
     delta_n = frec * n * (1-DefectDensity.M[l]/Na);
     #endif
     

     DefectDensity.M[l]
         = DefectDensity.M[l]+delta_n/volume;

if(DefectDensity.M[l]>Na) DefectDensity.M[l] = Na;  // Saturacion de la densidad de defectos

    }
   if(NNN%10000==0 && NNN!=0) WriteAmorphization();
   NNN++;
 }

 // Write amorphization histogram -----------------------------------------
void clase_Amorfizacion::WriteAmorphization()
  {
   ofstream Out;
   int m;
   char F[300];
   sprintf(F, "%s/Amorph1D.dat", DIRECTORY );
   Out.open(F);
   if(!Out) { cerr << "# I can't write in " << F << endl; exit(WRITE_ERROR); };
   Out << "# "<< DefectDensity.mPrimero << 
          " " << DefectDensity.mUltimo << 
          " " << DefectDensity.mAncho << endl;
   for(m=0; m<N_BINS; m++)
    {
     Out << m*DefectDensity.mAncho/10 << " "  // (nm)
     << DefectDensity.M[m]        << " "  // (defects/A^3)
     << DefectDensity.M[m]*100/Na << endl;// percentage
    }
   Out << "# Distance(nm) defects/A^3  percentage_of_amorphization\n"; 
   Out.close();
  }

 // SUMA Energia ----------------------------------------------------------
#ifdef ENERGIA_DEPOSITADA
void clase_Amorfizacion::Suma(double X, double Y, double Z, double Energy,
               char C[2], double Weigth)
  {
  // You need to use LevelOfSimulation=0 to simulate complete cascades and to proof
  // the total energy deposited in the crystal.
  // Total energy must to coincide with the product of the number of ions simulated 
  // by the initial energy of the ions.
  // 
  // R-> Remainder energy. Remainder energy of the projectile when we consider it is stopped
  // F-> Free flight energy loss
  // Q-> Local + Non local inelastic energy loss
  // B-> Nuclear transferred energy when an atom is going to be a projectile
  // I-> Nuclear transferred energy to an atom
   if(strcmp(C,"B")==0)  ETOTAL = ETOTAL + DisplacementEnergy*Weigth;
   else                  ETOTAL = ETOTAL + Energy;
   //OO << X << " " << Y << " " << Z << " " << Energy << " " << C << endl;
  };
#endif  

 // Destructor ------------------------------------------------------------
clase_Amorfizacion::~clase_Amorfizacion()
   {
#ifdef ENERGIA_DEPOSITADA
    cout << "####### Total energy = " << ETOTAL << endl;
#endif
   }

clase_Amorfizacion Amorfizacion;

}
// ------------------------------------------------------------------------
