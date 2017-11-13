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
// SIMUL.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "simul.h"

using std::ostream;
using std::endl;

namespace BCA {
ostream &operator<<( ostream &Out, SimulationDefinition A )
{
    int i_;
/*
    if(A.PROJS[0]==0)
     Out << _RED_ << _BOLD_ << "Implanted atom " 
         << _NORMAL_ << _BROWN_ << "= " 
         << _GREEN_  << A.AtomP << endl;
    else 
     {
      Out << _RED_ << _BOLD_ << "Implanted atoms ";
      for( i_=0; i_<N_ATOMS-1; i_++ ) Out << A.PROJS[i_] << ", ";
      Out << A.PROJS[i_] << endl;
     }  
    Out << _RED_ << _BOLD_ << "Energy         " 
        << _NORMAL_ << _BROWN_ << "= " 
        << _GREEN_  << A.ENERGY << " eV" << endl;
    if(A.Miller) 
      Out << _RED_ << _BOLD_ << "DIR            " 
          << _NORMAL_ << _BROWN_ << "= " 
          << _GREEN_ << A.DIR << "\t";
    else         
      Out << _RED_ << _BOLD_ << "Tha, Phi       " 
          << _NORMAL_ << _BROWN_ << "= " 
          << _GREEN_ << A.Tha << ", " << A.Phi << "\t";
    Out << _RED_ << _BOLD_ << "Divergency     " 
        << _NORMAL_ << _BROWN_ << "= " 
        << _GREEN_ << A.Divergency << endl;
    Out << _RED_ << _BOLD_ << "Implants       " 
        << _NORMAL_ << _BROWN_ << "= " 
        << _GREEN_ << A.NumberOfImplants << "\t\t";
    Out << _RED_ << _BOLD_ << "Dose           " 
        << _NORMAL_ << _BROWN_ << "= " 
        << _GREEN_ << A.Dose << "at/cm^2" << endl;
    Out << _RED_ << _BOLD_ << "DoseSplitting  ";
    if(A.DoseSplitting!=0)
    { 
     Out << _BOLD_ << _BROWN_ << "On" << endl;
     Out << _RED_ << _BOLD_ << "DoseDistribution " 
         << _NORMAL_ << _BROWN_ << "= " ;
     for(i_=0;i_<10;i_++) 
      if (A.DoseDistribution[i_]!=0) 
       Out << _GREEN_ << A.DoseDistribution[i_] << ",";
     Out << _GREEN_ << "0.0" << endl;
    }
    else Out << _BOLD_ << _BROWN_ << "Off" << _NORMAL_ << endl;
    Out << _RED_ << _BOLD_ << "Temperature    " 
        << _NORMAL_ << _BROWN_ << "= " 
        << _GREEN_  << A.Temperature << "\t\t";
    Out << _RED_ << _BOLD_ << "Rare Event     " 
        << _NORMAL_ << _BROWN_ << "= " 
        << _GREEN_ << A.SplitMag << endl;
    Out << _RED_ << _BOLD_ << "RS0            "
        << _NORMAL_ << _BROWN_ << "= " 
        << _GREEN_ ;
*/

  if(A.PROJS[0]==0)
     Out << "Implanted atom = " 
         <<  A.AtomP << endl;
    else 
     {
      Out << "Implanted atoms ";
      for( i_=0; i_<N_ATOMS-1; i_++ ) Out << A.PROJS[i_] << ", ";
      Out << A.PROJS[i_] << endl;
     }  
    Out << "Energy         = " 
        << A.ENERGY << " eV" << endl;
    if(A.Miller) 
      Out << "DIR            = " 
          << A.DIR << "\t";
    else         
      Out << "Tha, Phi       = " 
          << A.Tha << ", " << A.Phi << "\t";
    Out << "Divergency     = " 
        << A.Divergency << endl;
    Out << "Implants       = " 
        << A.NumberOfImplants << "\t\t";
    Out << "Dose           = " 
        << A.Dose << "at/cm^2" << endl;
    Out << "DoseSplitting  ";
    if(A.DoseSplitting!=0)
    { 
     Out << "On" << endl;
     Out << "DoseDistribution = " ;
     for(i_=0;i_<10;i_++) 
      if (A.DoseDistribution[i_]!=0) 
       Out << A.DoseDistribution[i_] << ",";
     Out << "0.0" << endl;
    }
    else Out << "Off" << endl;
    Out << "Temperature    = " 
        << A.Temperature << "\t\t";
    Out << "RS0            = " ;
    for(int k=0;k<MAXLAYER;k++) Out << A.RS0[k] << ",";
    Out << "0.0" << endl;
    Out << endl;
    return Out;
}
}
