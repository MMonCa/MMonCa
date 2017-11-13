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
// MAIN.CPP
//
// Version 2003.10.07
//
///////////////////////////////////////////////////////////////////////////

#include "main.h"

// Main program ////////////////////////////////////////////////////////////

using namespace BCA;
using std::cout;
using std::cerr;
using std::endl;

int main( int argc, char* argv[] )
{
 unsigned long int KK;
 int KK_Sim=1;
 NullVector(0,0,0);
 XVector(1,0,0);
 static MainLoop MAINLOOP;          

 cout << _CLS_    << endl 
      << _BOLD_   << _BLUE_ << "[IIS] " 
      << _NORMAL_ << _BLUE_ << "Ion Implant Simulator" 
      << _RED_ << " (Version 2003.10.07)" << _BOLD_  
      #ifdef DAMAGE1D      
      << " DAMAGE1D"  
      #endif
      #ifdef DAMAGE3D      
      << " DAMAGE3D"  
      #endif
      #ifdef DAMAGE3D_FASE1
      << " FASE1"  
      #endif
      #ifdef FAST_KP
      << " FAST_KP"  
      #endif
      << _NORMAL_ << endl;
 cout << _BLUE_ <<  "Developed at the Electronics Department of the "
      << "University of Valladolid, SPAIN" 
      << endl   << _RED_ << "(C) Dr. Jesus M. Hernandez-Mangas" << endl
      << _NORMAL_ << endl;
      
 MAINLOOP.Default();                    // Default values
 MAINLOOP.SetParameters( argc, argv );  // Read input file
 {
  MAINLOOP.Init( argc, argv );          // START

    #ifdef MODEL2
      FILE * O;
      char PP[300];
      sprintf( PP, "%s/MaxAmorf_vs_Dose.dat",DIRECTORY);
      O=fopen( PP,"w" );
      fclose(O);
    #endif

  LinkObject<SimulationDefinition>* Sim;
  Sim = SIMULACIONES.PointerToFirstData();
  do
  {

   MAINLOOP.InitSimulation( Sim->Data, KK_Sim );
  // 20100309 LuisJou MOD by J.H
  FILE *O_I,*O_V, *O_P, *O_PI, *O_BS;
  O_I  = fopen("Interstitials.txt","w");    
  O_V  = fopen("Vacancies.txt","w");    
  O_P  = fopen("Projectiles.txt","w");    
  O_PI = fopen("Projectiles_InitPos.txt","w");    
  O_BS = fopen("Backscattered.txt","w");    
  fprintf(O_I, "#             X              Y               Z                W          Energy Atom_index \n");
  fprintf(O_V, "#             X              Y               Z                W  Atom_index \n");
  fprintf(O_P, "#             X              Y               Z                W          Energy Atom_index \n");
  fprintf(O_PI,"#             X              Y               Z                W          Energy Atom_index \n");
  fprintf(O_BS,"#             X              Y               Z                W"
			   "           DirX           DirY            DirZ           Energy Atom_index \n");  
  fclose(O_I);
  fclose(O_V);
  fclose(O_P);
  fclose(O_PI);
  fclose(O_BS);
  //
   
   
#ifdef DAMAGE3D_FASE1   
   int TTTT;
   double PPPP= MAINLOOP.WinA*MAINLOOP.WinB*1e-16/ (MAINLOOP.NumberOfImplants/MAINLOOP.Dose);
   if(PPPP>=1.0) TTTT = (int) ceil(MAINLOOP.NumberOfImplants * PPPP) ;
   else          TTTT = MAINLOOP.NumberOfImplants;
   cerr << "####N_total= " << TTTT << endl;
#endif   

   do
   {
    MAINLOOP.DoIt();                    // DOIT cascade
    KK = MAINLOOP.DoItRareEvent(); 

  // 20100309 LuisJou MOD by J.H
  FILE *O_I,*O_V, *O_P;
  if(MAINLOOP.Extraccion & 0x01) 
  {
   O_I = fopen("Interstitials.txt","a+"); fprintf(O_I,"# Cascade %6ld\n\n",KK-1); fflush(O_I);
   fclose(O_I); 
  }
  if(MAINLOOP.Extraccion & 0x01) 
  {   
   O_V = fopen("Vacancies.txt","a+"); fprintf(O_V,"# Cascade %6ld\n\n",KK-1); fflush(O_V);  
   fclose(O_V);
  }
  
  if(MAINLOOP.Extraccion & 0x01) 
  {
   O_P = fopen("Projectiles.txt","a+"); fprintf(O_P,"# Cascade %6ld\n\n",KK-1); fflush(O_P);
   fclose(O_P);
  }
  
  if(MAINLOOP.Extraccion & 0x01) 
  {
   O_PI = fopen("Projectiles_InitPos.txt","a+"); fprintf(O_PI,"# Cascade %6ld\n\n",KK-1); fflush(O_PI);
   fclose(O_PI); 
  } 
   
  if(MAINLOOP.Extraccion & 0x01) 
  {
   O_BS = fopen("Backscattered.txt","a+"); fprintf(O_BS,"# Cascade %6ld\n\n",KK-1); fflush(O_BS);
   fclose(O_BS);
  }  

  if(MAINLOOP.Extraccion & 0x20)
   {
     FILE * O_MMONCA = fopen("Cascades.txt","a+"); fprintf(O_MMONCA, "# New cascade\n"); fflush(O_MMONCA);
     fclose(O_MMONCA);
   }
  //
    
    #ifdef MODEL2
    if(KK%1==0) 
     { 
      FILE * O;
      char PP[300];
      sprintf( PP, "%s/MaxAmorf_vs_Dose.dat",DIRECTORY);
      O=fopen( PP,"a" );
      #ifdef DAMAGE1D
      fprintf(O,"%8g %lf %lf %lf\n",
                (double) KK*MAINLOOP.Dose/MAINLOOP.NumberOfImplants, 
                (double) Amorfizacion.DefectDensity.Max()/Amorfizacion.Na*100,
                (double) Amorfizacion.DefectDensity.mAncho,
                (double) Amorfizacion.DefectDensity.Area() );
      #endif
      #ifdef DAMAGE3D
      fprintf(O,"%lf %lf\n",(double) KK, (double) Amorfizacion3D.DefectDensity.Max() );      
      #endif
      fclose(O);
     }
    #endif
    
   }
#ifndef DAMAGE3D_FASE1   
   while( KK <= MAINLOOP.NumberOfImplants ); 
#else
   while( KK <= TTTT ); 
#endif   
   MAINLOOP.Done(KK_Sim);               // END
   KK_Sim++;
   Sim = SIMULACIONES.PointerToNextData();
  } while(Sim!=NULL); 
 }
 MAINLOOP.JoinDoseSplitting();

// cout << "Fin!" << endl;
// scanf("%d",&KK);
}
///////////////////////////////////////////////////////////////////////////
