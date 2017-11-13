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
////////////////////////////////////////////////////////////////////////////
//
// SSF.CPP    Specific Screening Function
//
////////////////////////////////////////////////////////////////////////////
#include "ssf.h"

using std::cerr;
using std::ifstream;
using std::ofstream;
using std::endl;

namespace BCA {

void GetSSFName(int nZ1, int nZ2, char *Name )
{
 char C;
 
 switch( SpecificPotential[nZ1][nZ2] )
 { 
  case SPECIFIC     : C='S'; break; // Specific screening function
  case THOMAS_FERMI : C='T'; break; // Thomas-Fermi screening function
  case MOLIERE      : C='M'; break; // Moliere screening function
  case LENZ_JENSEN  : C='L'; break; // Lenz-Jensen screening function
  case BOHR         : C='B'; break; // Bohr screening function
  case SPECIFIC_WOPD: C='S'; break; // Specific ...
  default           : C='Z'; break; // ZBL screening funtion
 }
 sprintf(Name,"%cel%02d%02d.dat", C, nZ1, nZ2);
}

////////////////////////////////////////////////////////////////////////////

void SSF::Init( int Z1, int Z2, int Tipo )
{
    double X[N_SSF], Y[N_SSF], Y2[N_SSF];
    char F[200];
    char *TABLEDIR;
    double x,  S, DS;     
    int i;
    ifstream IN;
    
    _au = 0.8854 * 0.529 / (pow((double )Z1,0.23)+pow((double )Z2,0.23));
    switch(Tipo)
    {
     case SPECIFIC :    // Specific Screening function

        TABLEDIR = getenv( "TABLES3DDIR" );
        if( TABLEDIR == NULL ) strcpy( TABLEDIR, "./" );
        sprintf( F, "%sssf%02d%02d.dat", TABLEDIR, Z1, Z2 );
        IN.open( F );
        if( IN )
         {
          int i = 0,j;
          while( !IN.eof() )
          { 
           IN >> X[i];
           IN >> Y[i] >> Y2[i];
           i++;
          }
          int ND = i-1;
          for( x=0,i=0; x< D_SSF; x+= D_SSF/N_SSF,i++ )
           {
            int lo=1;
        int hi=ND;
        do
        {
         j=(hi+lo)/2;
         if( X[j]>x ) hi=j;
         else         lo=j;
        } while (hi-lo !=1);
        j=hi;
        double x0 = X[j-1];
        double x1 = X[j  ];
        double y0 = Y[j-1];
        double y1 = Y[j  ];
        y0 = (x-x0)*(y1-y0)/(x1-x0)+y0;
        _screen[i] = y0;
        y0 = Y2[j-1];
        y1 = Y2[j  ];
        y0 = (x-x0)*(y1-y0)/(x1-x0)+y0;
        _dscreen[i] = y0;
           }
         }
        else
         {
          cerr << "# Error: SSF::Init, I can't read the specific screening function " << F << endl;
          exit(READ_ERROR);
         }
        break;
     case SPECIFIC_WOPD :   // Specific Screening function WithOut Potential Derivative

        TABLEDIR = getenv( "TABLES3DDIR" );
        if( TABLEDIR == NULL ) strcpy( TABLEDIR, "./" );
        sprintf( F, "%sssf%02d%02d.wopd", TABLEDIR, Z1, Z2 );
        IN.open( F );
        if( IN )
         {
          int i = 0,j;
          while( !IN.eof() )
          { 
           IN >> X[i];
           IN >> Y[i];
           i++;
          }
          int ND = i-1;
          for( x=0,i=0; x< D_SSF; x+= D_SSF/N_SSF,i++ )
           {
            int lo=1;
        int hi=ND;
        do
        {
         j=(hi+lo)/2;
         if( X[j]>x ) hi=j;
         else         lo=j;
        } while (hi-lo !=1);
        j=hi;
        double x0 = X[j-1];
        double x1 = X[j  ];
        double y0 = Y[j-1];
        double y1 = Y[j  ];
        y0 = (x-x0)*(y1-y0)/(x1-x0)+y0;
        _screen[i] = y0;

           }
          double y0;
          y0=1;
          ofstream OUT;
          sprintf(F,"%sssf%02d%02d.w",TABLEDIR,Z1,Z2);
          OUT.open(F);
          for( x=0,i=0; x< D_SSF-D_SSF/N_SSF; x+= D_SSF/N_SSF, i++)
           {
        	  _dscreen[i]=(_screen[i]-y0)/(D_SSF/N_SSF);
            y0 = _screen[i];
            OUT << x << " " << _screen[i] << " " << _dscreen[i] << endl;
           }
         }
        else
         {
          cerr << "# Error, SSF::Init, I can't read the specific screening function " << F << endl;
          exit(READ_ERROR);
         }
        break;

    case THOMAS_FERMI :
        for( x=0,i=0; x< D_SSF; x+= D_SSF/N_SSF,i++ )
        {
            double X =x/_au;
            Thomas_Fermi_Screening( X, S, DS ); 
            _screen[i]  = S;
            //      DScreen[i] = DS/au;
        }   
        break;
    case MOLIERE :
        for( x=0,i=0; x< D_SSF; x+= D_SSF/N_SSF,i++ )
        {
            double X =x/_au;
            Moliere_Screening( X, S, DS );  
            _screen[i]  = S;
            _dscreen[i] = DS/_au;
        }   
        break;
    case LENZ_JENSEN :
        for( x=0,i=0; x< D_SSF; x+= D_SSF/N_SSF,i++ )
        {
            double X =x/_au;
            Lenz_Jensen_Screening( X, S, DS );  
            _screen[i]  = S;
            _dscreen[i] = DS/_au;
        }   
        break;
    case BOHR :
        for( x=0,i=0; x< D_SSF; x+= D_SSF/N_SSF,i++ )
        {
            double X =x/_au;
            Bohr_Screening( X, S, DS ); 
            _screen[i]  = S;
            _dscreen[i] = DS/_au;
        }   
        break;
    default :   // ZBL Screening function
        for( x=0,i=0; x< D_SSF; x+= D_SSF/N_SSF,i++ )
        {
            double X =x/_au;
            ZBL_Screening( X, S, DS );
            _screen[i]  = S;
            _dscreen[i] = DS/_au;
       }
    }    
}

////////////////////////////////////////////////////////////////////////////

void SSF::ZBL_Screening( double X, double &S, double &DS )
    {
        double A[4];
        A[0] = -3.2; A[1] = -0.9423; A[2] = -0.4029; A[3] = -0.2016;
        double B[4];
        B[0] =  0.18179; B[1] = 0.50986; B[2] = 0.28018; B[3] = 0.02817;
    double EXP;

    EXP = exp( A[0] * X );
    S   = B[0] * EXP;
    DS  = A[0]*B[0] * EXP;
    EXP = exp( A[1] * X );
    S      += B[1] * EXP;
    DS     += A[1]*B[1] * EXP;
    EXP = exp( A[2] * X );
    S      += B[2] * EXP;
    DS     += A[2]*B[2] * EXP;
    EXP = exp( A[3] * X );
    S      += B[3] * EXP;
    DS     += A[3]*B[3] * EXP;
    }
    
void SSF::Thomas_Fermi_Screening( double X, double &S, double &DS )
    {
        double l  = 0.8034;
    double a  = pow(12.0,0.6666666);
        S       = pow( 1+pow(X/a,l) ,-3/l);
//  DS      = -3/l*pow(1+pow(X/a,l),-(3+l)/l)*l/pow(a,l)*pow(X,l-1); //ERROR
    }
    
void SSF::Moliere_Screening(double X, double &S, double &DS )
    {
        double A[3];
        A[0] = -0.30; A[1] = -1.20; A[2] = -6.0;
        double B[3];
        B[0] =  0.35; B[1] =  0.55; B[2] =  0.1;
        
    double EXP;

    EXP = exp( A[0] * X );
    S   = B[0] * EXP;
    DS  = A[0]*B[0] * EXP;
    EXP = exp( A[1] * X );
    S      += B[1] * EXP;
    DS     += A[1]*B[1] * EXP;
    EXP = exp( A[2] * X );
    S      += B[2] * EXP;
    DS     += A[2]*B[2] * EXP;
    }
    
void SSF::Lenz_Jensen_Screening(double X, double &S, double &DS )
    {
        double A[3];
        A[0] = -1.0380; A[1] = -0.3876; A[2] = -0.2060;
        double B[3];
        B[0] =  0.7466; B[1] =  0.2433; B[2] =  0.01018;
    double EXP;

    EXP = exp( A[0] * X );
    S   = B[0] * EXP;
    DS  = A[0]*B[0] * EXP;
    EXP = exp( A[1] * X );
    S      += B[1] * EXP;
    DS     += A[1]*B[1] * EXP;
    EXP = exp( A[2] * X );
    S      += B[2] * EXP;
    DS     += A[2]*B[2] * EXP;
    }
    
void SSF::Bohr_Screening(double X, double &S, double &DS )
    {
        S       =  exp( -X );
    DS      = -exp( -X );
    }

double SSF::Screening( double R, double &S, double &DS )
    {
    // Input    R : Distance (A)
    // Output   S : Screening (adimensional)
    //      DS: Derivative of screening (adimensional)

    double F1=0, F2=0;

    int i  = (int) (R / (D_SSF/N_SSF));
    if ( i < 0 )
    {
    	std::cout << "i < 0\n";
    	//exit(1);
    }
    if( i >= N_SSF - 2 )
     {
        S = 0;
        DS = 0;
     }
    else
     {
        F1 = _dscreen[i];
        F2 = _dscreen[i+1];

        double K = ( R/(D_SSF/N_SSF) - i );
        DS = F1 + (F2-F1)*K;

        F1 = _screen[i];
        F2 = _screen[i+1];
        S  = F1 + (F2-F1)*K;
     }
    return S;
    }

double SSF::ScreeningOnly( double R, double &S )
    {
    // Input    R : Distance (A)
    // Output   S : Screening (adimensional)
    int i;
    double F1, F2;

    i  = (int) (R / (D_SSF/N_SSF));
    if( i >= N_SSF - 2 )
    {
        S = 0;
        return 0;
    }
    else
     {
        F1 = _screen[i];
        F2 = _screen[i+1];
        S = F1 + (F2-F1)*( R/(D_SSF/N_SSF) - i );
        return S;
     }
    }

}
////////////////////////////////////////////////////////////////////////////
