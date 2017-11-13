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
// POTEN.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "poten.h"
namespace BCA {
ZBL::ZBL()  // Constructor
{           // Initialize table
        // Coefs for ZBL potential
        double _EXPAR[4] = { -3.2    , -0.9423,  -0.4029, -0.2016  };
        double  _VPAR[4] = {  0.18179,  0.50986,  0.28018, 0.02817 };

        register int  j;

        // Initialize coefs.
        for( j=0; j<4; j++ )
         {
            EXPAR[ j ] =  _EXPAR[ j ];
             VPAR[ j ] =   _VPAR[ j ];
            DVPAR[ j ] = -_EXPAR[ j ] * _VPAR[ j ];
         };

        // Initialize cache
 #ifdef CACHED_ZBL
        double X;

        X  = 0;
        DX = (double ) NUMMAX/MATR;

        for( j=0; j<MATR; j++ )
         {
           double EXP, VR, DVR;

           EXP  = exp( EXPAR[0] * X );
           VR   =  VPAR[0] * EXP;
           DVR  = DVPAR[0] * EXP;
           EXP  = exp( EXPAR[1] * X );
           VR  +=  VPAR[1] * EXP;
           DVR += DVPAR[1] * EXP;
           EXP  = exp( EXPAR[2] * X );
           VR  +=  VPAR[2] * EXP;
           DVR += DVPAR[2] * EXP;
           EXP  = exp( EXPAR[3] * X );
           VR  +=  VPAR[3] * EXP;
           DVR += DVPAR[3] * EXP;

           TABLEV [j] = VR;
           TABLEDV[j] = DVR;

           X = X + DX;
         };
 #endif
}

#ifdef CACHED_ZBL
void ZBL::Potential( double Xa, double &VR, double &DVR )
     {
        double XX;
        int  X1;
        double EXP;

        if( Xa < NUMMAX - DX )  // With interpolation
        {

          XX = Xa *MATR/NUMMAX;
          X1 = (int) XX;

          XX = XX-X1;

          VR = TABLEV [X1] + ( TABLEV [X1+1] - TABLEV [X1] )* XX;
         DVR = TABLEDV[X1] + ( TABLEDV[X1+1] - TABLEDV[X1] )* XX;
        }
        else            // Without interpolation
        {
         EXP  = exp( EXPAR[0] * Xa );
         VR   =  VPAR[0] * EXP;
         DVR  = DVPAR[0] * EXP;
         EXP  = exp( EXPAR[1] * Xa );
         VR  +=  VPAR[1] * EXP;
         DVR += DVPAR[1] * EXP;
         EXP  = exp( EXPAR[2] * Xa );
         VR  +=  VPAR[2] * EXP;
         DVR += DVPAR[2] * EXP;
         EXP  = exp( EXPAR[3] * Xa );
         VR  +=  VPAR[3] * EXP;
         DVR += DVPAR[3] * EXP;
        };
};

void ZBL::PotentialOnly( double Xa, double &VR )
{
        double XX;
        int  X1;
        double EXP;

        if( Xa < NUMMAX - DX )  // With interpolation
        {

          XX = Xa *MATR/NUMMAX;
          X1 = (int) XX;

          XX = XX-X1;

          VR = TABLEV [X1] + ( TABLEV [X1+1] - TABLEV [X1] )* XX;
        }
        else            // Without interpolation
        {
         EXP  = exp( EXPAR[0] * Xa );
         VR   =  VPAR[0] * EXP;
         EXP  = exp( EXPAR[1] * Xa );
         VR  +=  VPAR[1] * EXP;
         EXP  = exp( EXPAR[2] * Xa );
         VR  +=  VPAR[2] * EXP;
         EXP  = exp( EXPAR[3] * Xa );
         VR  +=  VPAR[3] * EXP;
        };
};

#else   // Uncached function
void ZBL::Potential( double Xa, double &VR, double &DVR )
     {
        double EXP;

        EXP  = exp( EXPAR[0] * Xa );
        VR   =  VPAR[0] * EXP;
        DVR  = DVPAR[0] * EXP;
        EXP  = exp( EXPAR[1] * Xa );
        VR  +=  VPAR[1] * EXP;
        DVR += DVPAR[1] * EXP;
        EXP  = exp( EXPAR[2] * Xa );
        VR  +=  VPAR[2] * EXP;
        DVR += DVPAR[2] * EXP;
        EXP  = exp( EXPAR[3] * Xa );
        VR  +=  VPAR[3] * EXP;
        DVR += DVPAR[3] * EXP;
}

void ZBL::PotentialOnly( double Xa, double &VR )
{
        double EXP;

        EXP  = exp( EXPAR[0] * Xa );
        VR   =  VPAR[0] * EXP;
        EXP  = exp( EXPAR[1] * Xa );
        VR  +=  VPAR[1] * EXP;
        EXP  = exp( EXPAR[2] * Xa );
        VR  +=  VPAR[2] * EXP;
        EXP  = exp( EXPAR[3] * Xa );
        VR  +=  VPAR[3] * EXP;
}
}
#endif

///////////////////////////////////////////////////////////////////////////
