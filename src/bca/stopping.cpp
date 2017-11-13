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
// STOPPING.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "stopping.h"

using namespace std;
namespace BCA {
// External variables ----------------
extern int LayerNumber;

// Global variables ------------------
STP STPO;
double rs0;
double RS0_LAYERS[MAXLAYER];
int Frenado=2;
int Ionizacion=1;

struct StoppingStruct ESA[ ] = {
                { "None"    ,"No stopping",             NoStopping  },
                { "ZBL"     ,"BK stopping",             ZBLStopping },
                { "Our"     ,"Our model",               OurStopping },
                { "?"       ,"This help",               NoStopping  }
                   };

struct IonizationStruct IS[ ] = {
                { "ZBL"     ,"ZBL ionization", ZBL_I },
                { "BK"      ,"BK  ionization",  BK_I },
                { "CGJ"     ,"CGJ ionization", CGJ_I },
                { "MP"      ,"MP  ionization",  MP_I }
                };

// Ionization curves //////////////////////////////////////////////////////

double ZBL_I( double yr )    // => ZIEGLER, BIERSAK, LITMARK
{
 double q=1.-exp(0.803*pow(yr,0.3)-1.3167*pow(yr,0.6)-0.38157*yr-0.008983*yr*yr);
 if (q<=0.0) return 0.0;
 else        return q;
}

double CGJ_I( double yr )    // => CAI, GROENBECH-JENSEN
{
 double q=1.-exp(-0.95*(yr-0.07));
 if (q<=0.0) return 0.0;
 else        return q;
}

double BK_I( double yr )     // => BRANDT, KITAGAWA
{
 double q=1.-exp(-0.92*yr);
 if (q<=0.0) return 0.0;
 else        return q;
}

double MP_I( double yr )     // => MATHAR, POSSELT
{
 double q;
 if(yr<0.6) q = 1.1111111*yr*yr; 
 //q=1-exp(-0.95*(yr-0.0700));
 else       q=1-exp(-2.00*(yr-0.3483));
 if (q<=0.0) return 0.0;
 else        return q;
}

// Stopping functions /////////////////////////////////////////////////////

double  NoStopping ( double v, double z, double rs )
 {
  return( 0 );
 }

double  ZBLStopping(double v,double z, double rs)
{
    double vf,f,fr,fb,s0,yr,vr,vmb,wp,q,a,gi,x,stopping;

    vf= 1.919158292678/rs;
    x=0.1658591099016*rs;
    f=2./3./M_PI*(log((1.+2.*x/3.)/x)-(1.-x/3.)/(1.+2.*x/3.))*vf
                            /(1.-x/3.)/(1.-x/3.);
    fr=f*(1.644*exp(-((rs-2.197)*(rs-2.197))/15.8));
    wp=sqrt(3./rs/rs/rs);
    fb=vf/v* wp*wp/v/v*log(2*v*v/wp);
    vmb=sqrt(wp/2*M_E);
    if (v>=vmb)  // ZONE II // Get the lower curve
       if (fr>=fb)  f = fb;
       else     f = fr;
    else         // ZONE I  // Get the linear aproximation
       f = fr;
    f=f*51.42;
    s0=v/vf*f;  // (eV/A)
    if (z<2.) return(s0);
    if (v<=vf)
        vr=3./4.*vf*(1.+2./3.*v*v/vf/vf-1./15.*v*v*v*v/vf/vf/vf/vf);
    else
        vr=v*(1.+1./5.*vf*vf/v/v);
    yr=vr/pow(z,0.6666666);     //  vb = 1 in atomic units

    q = IS[Ionizacion].Ionization( yr );

    a=(0.4801*(pow((1.-q),0.666666)))/((pow(z,0.333333))*(1.-0.143*(1.0-q)));
    gi=q+0.5*(1.0-q)* log(1.0+(4.0*a/rs)*(4.0*a/rs));
    stopping=gi*z*gi*z*s0;
    return(stopping);
}

double  OurStopping(double v,double z, double rs)
{
    double vf,f,fr,fb,s0,yr,vr,vmb,wp,q,a,gi,x,stopping;

    vf= 1.919158292678/rs;
    x=0.1658591099016*rs;
    f=2./3./M_PI*(log((1.+2.*x/3.)/x)-(1.-x/3.)/(1.+2.*x/3.))*vf
                            /(1.-x/3.)/(1.-x/3.);
    fr=f*(1.644*exp(-((rs-2.197)*(rs-2.197))/15.8));

    wp=sqrt(3./rs/rs/rs);
    fb=vf/v* wp*wp/v/v*log(2*v*v/wp);
    vmb=sqrt(wp/2*M_E);
    if (v>=vmb)  // ZONE II // Get the lower curve
       if (fr>=fb)  f = fb;
       else     f = fr;
    else         // ZONE I  // Get the linear aproximation
       f = fr;
    f=f*51.42;
    s0=v/vf*f;  // (eV/A)
    if (z<2.) return(s0);
    
    rs= rs0; //(a.u)

    vf= 1.919158292678/rs;
    
    if (v<=vf)
        vr=3./4.*vf*(1.+2./3.*v*v/vf/vf-1./15.*v*v*v*v/vf/vf/vf/vf);
    else
        vr=v*(1.+1./5.*vf*vf/v/v);
    yr=vr/pow(z,0.6666666);     //  vb = 1 in atomic units

    q = IS[Ionizacion].Ionization( yr );

    a=(0.4801*(pow((1.-q),0.666666)))/((pow(z,0.333333))*(1.-0.143*(1.0-q)));
    gi=q+0.5*(1.0-q)* log(1.0+(4.0*a/rs)*(4.0*a/rs));
    stopping=gi*z*gi*z*s0;
    return(stopping);
}

///////////////////////////////////////////////////////////////////////////
//
//  STP CLASS
//

STP::STP()
{
  int i,j;
  for(i=0;i<104;i++)
     for(j=0;j<5;j++) TABLA[i][j] = NULL;
}

void STP::Clear()
{
  int i,j;
  for(i=0;i<104;i++)
   for(j=0;j<5;j++)
    if(TABLA[i][j]!=NULL)
    {
     free(TABLA[i][j]);
     TABLA[i][j]=NULL;
    }
}
         
double STP::Stopping( double V, double dZ, double rS )
{
 double *NT;
 double v,rs;
 int i,j;
 int Z = (int) dZ;

 if( TABLA[Z][LayerNumber]==NULL )
  {
   char FF[300];
   ifstream IF;
   ofstream OF;

   NT = (double*) malloc( NTABLEX*NTABLEY*sizeof(double) );
   if( NT== NULL ){ cerr << "#Error: STP::Stopping, No memory\n"; exit(MEMORY_ERROR); };

   char * TABLEDIR, Directory[255];   
   TABLEDIR = getenv("TABLES3DDIR");         
   if( TABLEDIR != NULL ) strcpy( Directory, TABLEDIR );
   else                   strcpy( Directory, TABLEDIRDEFAULT );
   sprintf(FF,"%sInelNonLocal%02d_%02d_%02d_%04.2f.dat", Directory,
                        Z, Frenado, Ionizacion, (float)rs0);
   #ifdef WINDOWS
   IF.open(FF,ios::binary|ios::nocreate);
   #else
   IF.open(FF);
   #endif
   if(IF)
    {
     cerr << "# Reading " << FF ;
#ifdef GCC3 
     IF.read( (char*) NT,(int) NTABLEX*NTABLEY*sizeof(double) );
#else
     IF.read( (unsigned char*) NT, NTABLEX*NTABLEY*sizeof(double) );
#endif 
    IF.close();
     TABLA[Z][LayerNumber] = NT;
     cerr << _OK_ << endl;
    }
   else
    {
     cerr << "# Generating non-local inelastic stopping table Z = " << Z;
     // Genero la tabla
     for(i=0,rs=0.001; i<NTABLEX; rs+=DLIMX/NTABLEX,i++)
      for(j=0, v =0.001; j<NTABLEY;  v+=DLIMY/NTABLEY,j++)
       {
        NT[i*NTABLEY+j] = ESA[Frenado].Stopping(v, dZ, rs);
       }
     TABLA[Z][LayerNumber] = NT;
     cerr << _OK_ << endl;
     #ifdef WINDOWS
     OF.open(FF,ios::binary);
     #else
     OF.open(FF);
     #endif
     if(!OF) {
      cerr << "# Error, STP::Stopping, Can't write in " << FF << endl;
      exit(WRITE_ERROR);
     };
#ifdef GCC3
     OF.write( (char*) NT, (int) NTABLEX*NTABLEY*sizeof(double) );
#else
     OF.write( (unsigned char*) NT, NTABLEX*NTABLEY*sizeof(double) );
#endif     
     OF.close();
    }
  }
 // Devuelvo el valor
 double F1,F2,F3,F4;
      
 i = (int) (rS/(DLIMX/NTABLEX));
 j = (int) ( V/(DLIMY/NTABLEY));

 if(i>NTABLEX-2||j>NTABLEY-2)
  {
   F1 = ESA[Frenado].Stopping(V,dZ,rS);
   return F1;
  }
    
 F1 = TABLA[Z][LayerNumber][ i   *NTABLEY + j   ];
 F2 = TABLA[Z][LayerNumber][ i   *NTABLEY +(j+1)];
 F3 = TABLA[Z][LayerNumber][(i+1)*NTABLEY + j   ];
 F4 = TABLA[Z][LayerNumber][(i+1)*NTABLEY +(j+1)];
    
 F1 = F1 +(F2-F1)*(  V*NTABLEY/DLIMY - j );    
 F3 = F3 +(F4-F3)*(  V*NTABLEY/DLIMY - j );

 F1 = F1 +(F3-F1)*( rS*NTABLEX/DLIMX - i );
 return F1;
}
}
