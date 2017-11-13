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
///////////////////////////////////////////////////////////////////////////f
//
//  INEL.CPP
//  Inelastic stopping due to electron-electron interaction
//
///////////////////////////////////////////////////////////////////////////
#include "inel.h"

using std::cerr;
using std::ofstream;
using std::ifstream;
using std::endl;

namespace BCA {

double InelT::FFF(  double R, double Z1, double Z2 )
  {  // R in (A)
   double *NT;
   double d, F, Z;
   int i;
   if ( Z2>Z1 ) { Z=Z1;Z1=Z2;Z2=Z; } 
   
   if( TABLA[(int)Z1][(int)Z2]==NULL )
    { 
     Ssf.Init((int)Z1,(int)Z2,ZBL_SCREENING);
     au = 0.8854*0.529/(pow(Z1,0.23)+pow(Z2,0.23));
     NT = (double *)malloc(NTABLE*sizeof(double));
     if( NT==NULL )
      { cerr << "# Error: Inel::FFF, No memory\n"; exit(MEMORY_ERROR); }
    
     ifstream IF;
     ofstream OF;
     char FF[300];

     char * TABLEDIR, Directory[255];   
     TABLEDIR = getenv("TABLES3DDIR");         
     if( TABLEDIR != NULL ) strcpy( Directory, TABLEDIR );
     else                   strcpy( Directory, TABLEDIRDEFAULT );

     sprintf(FF,"%sInelLocal_%02d_%02d.dat",Directory,(int)Z2,(int)Z1); 
     #ifdef WINDOWS
     IF.open(FF,ios::binary|ios::nocreate);
     #else
     IF.open(FF);
     #endif
     if(IF)
      {
       cerr << "# Reading " << FF;
#ifdef GCC3
       IF.read((char*)NT,NTABLE*sizeof(double));
#else       
       IF.read((unsigned char*)NT,NTABLE*sizeof(double));
#endif
       IF.close();
       TABLA[(int)Z1][(int)Z2] = NT; 
       cerr  << _OK_ << endl;       
      }
     else
      { 
       cerr << "# Generating     local inelastic stopping table, Z1 = "
            << Z1 << ", Z2 = " << Z2 << " : ";
       double alfa = 1/( 1 + pow( Z2/Z1, 1./6.)  );
       double a    = 0.4685025144;  // pow(9*M_PI*M_PI/128.,1./3.)*a0;
       for( i=0,d = 0.001; d <= DLIM; d += DLIM/NTABLE,i++ )
        {
          NT[i] = 7.549079e-15 * (Z1*Z1* I( pow(Z1,1./3.)*alfa*d/a   )
                            +Z2*Z2*I( pow(Z2,1./3.)*(1-alfa)*d/a  )  );
        }     
       TABLA[(int)Z1][(int)Z2] = NT; 
       cerr << i  << _OK_ << endl;
       #ifdef WINDOWS
       OF.open(FF,ios::binary);
       #else
       OF.open(FF);
       #endif
#ifdef GCC3
       if(!OF ||
          !OF.write((char*)NT,NTABLE*sizeof(double)))
#else
       if(!OF ||
          !OF.write((unsigned char*)NT,NTABLE*sizeof(double)))
#endif
       {
        cerr << "# Error: Inel::FFF, Can't write in " << FF << endl; 
        exit(WRITE_ERROR);
       };
       OF.close();
      }
    }

     double F1, F2;
     i = (int) (R /(DLIM/NTABLE));   
     if(i>NTABLE-1)
      { 
        double alfa = 1/( 1 + pow( Z2/Z1, 1./6.)  );
        double a    = 0.4685025144;
        F = 7.549079e-15 * (Z1*Z1* I( pow(Z1,1./3.)*alfa*R/a   )
                        +Z2*Z2*I( pow(Z2,1./3.)*(1-alfa)*R/a  )  );
        return F;   
      };
     F1 = TABLA[(int)Z1][(int)Z2][i];
     F2 = TABLA[(int)Z1][(int)Z2][i+1];    
     F = F1 + (F2-F1)*( R*NTABLE/DLIM - i);
     return F;
  }

InelT IT;       // Declare a instance

// ------------------------------------------------------------------------
// Integration in a trajectory

double InelLosses( double Z1, double M1, Vector PAbs, Vector Dir,
           double ENERGY, double E1, double AllowedDistance,
           double Z2, Vector TAbs, Vector TAbs1, Vector TDir1,
           Vector DirFinal )
  {
   Vector DR, TR, Delta, DeltaT;
   double DeltaR;
   double V, Fr, Distance, q, q2, DELTA;

   DR      = PAbs;                          // (angstroms = A)
   TR      = TAbs;                          // (A)
   DELTA   = AllowedDistance / NPASOS_ ;    // (A)
   Delta   = Dir*DELTA ;                    // (A)
   DeltaT  = 1.0/NPASOS_*( TAbs1-TAbs );    // (A)
   Distance = 0;
   q = 0.0;
   double CTE1 = sqrt( CTE/M1 );
   double CTE2 = (E1-ENERGY)/AllowedDistance;

// INCOMING trajectory -----------------------------------
   while( Distance < AllowedDistance )
   { 
    V = CTE1*sqrt( ENERGY + CTE2*Distance );
    DeltaR = (double ) (DR - TR);
    Fr = IT.FFF( DeltaR, Z1, Z2 );
    DR = DR + Delta;
    TR = TR + DeltaT;

    double VvB = V*bohr_velocity;
    double Einel1 =  Fr * VvB;
    double Einel2 =  Fr * Vcritica / VvB;
// double F_v =
//  (2*exp(-(V*V/Vcritica/Vcritica)))/(1+exp(-2*(V*V/Vcritica/Vcritica)));
    double expA = exp( -V*V/(Vcritica*Vcritica) );
    double F_v = 2*expA/(1+expA*expA);
// double Einel = F_v*Einel1 + (1-F_v)*Einel2;
    double Einel = F_v*(Einel1-Einel2) + Einel2;
    q += Einel;
    Distance += DELTA;
   };

// --------------------------------------   

   Delta = DirFinal * DELTA;
   DeltaT = Unitary(TDir1)*(double )DeltaT;
   Distance = 0;
   q2 = 0.0; 
   V  = CTE1*sqrt( ENERGY );

   double cV1 = V * bohr_velocity;
   double cV2 = Vcritica / cV1;
   double expA = exp( -V*V/(Vcritica*Vcritica) );
   double F_v = 2*expA/(1+expA*expA);

// OUTCOMING trajectory ------------------------------------
   while( Distance <  AllowedDistance )
   { 
    DeltaR = (double ) (DR - TR);
    Fr = IT.FFF( DeltaR, Z1, Z2 );
    DR = DR + Delta;
    TR = TR + DeltaT;
    double Einel = F_v*Fr*(cV1-cV2) + Fr*cV2;
    q2 += Einel;
    Distance += DELTA;
   };
   return (q+q2)*DELTA*1e-10/electron_charge;   // Put in correct units
  }
// ------------------------------------------------------------------------
}
