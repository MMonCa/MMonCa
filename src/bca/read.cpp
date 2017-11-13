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
// READ.CPP
//
///////////////////////////////////////////////////////////////////////////
#include "read.h"

using std::endl;
using std::cerr;
using std::ostream;
using std::istream;
using std::ofstream;
using std::ifstream;

namespace BCA {
// Output stream overloaded
ostream &operator<<( ostream &Out, Variable A)
 {
   int i,j,k;
          int  *cI;
          double *cR;
          char *cC;
        Vector *cV;
          long *cL;
  List<AtomDefinition> *cA;
  List<XtalDefinition> *cX;
 List<LayerDefinition> *cY;
 List<SimulationDefinition> *cS;
    short unsigned int *cIA;
//  Out << _BLUE_ << _BOLD_;
  Out.width(25);
//  Out << A.Name << _NORMAL_ << _BROWN_ << " = " << _GREEN_;
  Out << A.Name << " = ";
  Out.width(10);
  Out.setf(std::ios::showpoint);
  Out.precision(4);

  switch (A.Type)
  {
   case INTEGER: cI = (   int*) A.Address; Out << *cI << endl; break;
   case REAL   : cR = (  double*) A.Address; Out << *cR << endl; break;
   case CHAR   : cC = (  char*) A.Address; Out << *cC << endl; break;
   case TEXT_   : cC = (  char*) A.Address;
            //do { Out << *cC; *cC++;} while ( *cC!= '\0' );
            Out << (const char *)A.Address << endl;
            break;
   case VECTOR : cV = (Vector*) A.Address; Out << *cV << endl; break;
   case LONG   : cL = (  long*) A.Address; Out << *cL << endl; break;
   case ATOMDEFINITION:
     cA = (List<AtomDefinition>*) A.Address; Out << *cA << endl; break;
   case XTALDEFINITION:
     cX = (List<XtalDefinition>*) A.Address; Out << *cX << endl; break;
   case XTALDEFINITION2:
     cX = (List<XtalDefinition>*) A.Address; Out << *cX << endl; break;
   case LAYERDEFINITION:
    cY = (List<LayerDefinition>*) A.Address; Out << *cY << endl; break;
   case SIMULATIONDEFINITION:
    cS = (List<SimulationDefinition>*) A.Address; Out << *cS << endl; break;
   case ARRAY  : 
        cIA = (short unsigned int *) A.Address; 
    Out << endl;
    for(i=0;i<100;i++)
     {
     for(j=0;j<100;j++)
       if( *(cIA+i*100+j)!=0 )  Out << " Z1 = " << i 
                                   << ", Z2 = " << j 
                       << ", Type = " << *(cIA+i*100+j)<< endl;
     } 
    break;
   case MATRIX:
   case MATRIX2 :
        cR = (  double*) A.Address;
        for(k=0;k<MAXLAYER;k++) Out << *cR++ << ",";
        Out << endl;
        break;
   case INTMATRIX :
        cI = (   int*) A.Address;
        for(k=0;k<N_ATOMS;k++) Out << *cI++ << ",";
        Out << endl;
        break;
   case BOOLEAN :
        cI = (   int*) A.Address; 
        if( *cI!=0 ) Out << "On" << endl;
        else         Out << "Off" <<endl;
        break;
  };
//  Out << _NORMAL_;
  return Out;
 }

///////////////////////////////////////////////////////////////////////////

// Read parameters from a file ( return -1 if error )
int Input::ReadParameters( int TIP )
{
    // Local variables
    int      i,j,k,VV;        // Counters
    int      NVars = 0;     // Number of variables readed

          char CAR[255];
          double Value;   // Value readed
        Vector V;
    AtomDefinition A;
    XtalDefinition X;
   LayerDefinition* Y;
   SimulationDefinition* S=NULL;

    LinkObject<Variable> *Temp; // Temporary variable
            int  *cI;
            double *cR;
            char *cC;
          Vector *cV;
            long *cL;
    List<AtomDefinition> *cA;
    List<XtalDefinition> *cX;
   List<LayerDefinition> *cY;
   List<SimulationDefinition> *cS;
     short unsigned int  *cAI;
            
    std::string Name = "String readed";    // String readed

    // Open the stream
    ifstream IN( FileName.c_str() );
    if (!IN )
     {
      cerr << "# Error: Input::ReadParameters, Can't open " << FileName << endl;
      exit(READ_ERROR);
     };

    // Read the parameters
    while (!IN.eof() && Name[0]!=0)
    {
      IN >> Name;
      j = 0;
      if( !CaseSensitive )
       while ( Name[j] != 0 )
    { Name[j] = (char) _toupper( Name[j] ); j++; };

      Temp = PointerToFirstData();
      for( i=0; i<NumberOfElements(); i++)
       {
    j = 0;
    if( !CaseSensitive )
     while ( Temp->Data.Name[j] != 0 )
      { Temp->Data.Name[j] = (char) _toupper( Temp->Data.Name[j] ); j++; };

    if( Temp->Data.Name == Name)
    {
      switch (Temp->Data.Type)
      {
       case INTEGER : IN >> Value;
               cI = (int* ) Temp->Data.Address;
              *cI = (int  ) Value;
              break;
       case REAL    : IN >> Value;
               cR = (double*) Temp->Data.Address;
              *cR =     Value;
              break;
       case CHAR    : IN >> CAR;
               strcpy((char*)Temp->Data.Address,CAR);
              break;
       case TEXT_    : cC = (char*) Temp->Data.Address;
              strcpy( cC, "" );
              IN >> CAR;
                  while( CAR[0] !=']' )
                  {
               strcat( cC, CAR );
               strcat( cC, " " );
               if( IN.eof() )
                { 
                  cerr << "# Error: Input::ReadParameters, Incomplete string in " << FileName << "\n";
                  exit(INCOMPLETE_STRING_ERROR);
                }
               IN >> CAR;
              };
                      break;          
       case VECTOR  : IN >> V;
               cV = (Vector*) Temp->Data.Address;
              *cV =     V;
              break;
       case LONG    : IN >> Value;
               cL = (long*) Temp->Data.Address;
              *cL = (long ) Value;
              break;
       case ATOMDEFINITION :
              IN >> A;
               cA = (List<AtomDefinition>*) Temp->Data.Address;
              (cA)->Add( (AtomDefinition) A);
              break;
       case XTALDEFINITION :
              int tmp;
              IN >> X.Type >> tmp >> X.Origin >> X.BindingEnergy;
              X.Center = (CenteredType) tmp;
              X.Type2 = X.Type; 
              X.X = 1.0;
              X.BindingEnergy2 = X.BindingEnergy;
               cX = (List<XtalDefinition>*) Temp->Data.Address;
              (cX)->Add( (XtalDefinition) X);
              break;
       case XTALDEFINITION2 :
              IN >> X;
               cX = (List<XtalDefinition>*) Temp->Data.Address;
              (cX)->Add( (XtalDefinition) X);
              break;
       case LAYERDEFINITION :
              // Without parameters
                Y = (LayerDefinition*) Temp->Data.AuxiliaryAddress;
               cY = (List<LayerDefinition> *) (Temp->Data.Address);
               cY->Add( (LayerDefinition)*Y );
                Y->Default(); // Initialize XTal list
              break;
       case SIMULATIONDEFINITION :
              // Without parameters
                S = (SimulationDefinition*) Temp->Data.AuxiliaryAddress;
               cS = (List<SimulationDefinition>*) Temp->Data.Address;
              (cS)->Add( (SimulationDefinition) *S);
              break;
       case ARRAY : int I,J,K;
                       IN >> I >> J >> K;
               cAI= (short unsigned int*) Temp->Data.Address;
               *(cAI+I*100+J) = K;
               break;
       case MATRIX  :
               VV=0;
               cR = (double*) Temp->Data.Address;
               for(k=0;k<=MAXLAYER;k++) 
               {
                IN >> Value;
                if(Value!=0) *cR++ = Value;
                else 
                 {
                  VV = k;
                  Value = *(cR-1);
                  for(k=k;k<MAXLAYER;k++) *cR++ = Value;
                  break;
                 }
               }
               if(k==MAXLAYER+1&&VV==0) 
                {
                    cerr << "# Error: Input::ReadParameters, " << Temp->Data.Name
                         << " token must end with 0.0" << endl;
                    cerr << "MAXLAYER = " << MAXLAYER << endl;                         
                    exit(NOT_END_WITH_ZERO_ERROR);
                };
              break;
       case MATRIX2  :
               VV=0;
               cR = (double*) Temp->Data.Address;
               for(k=0;k<=MAXLAYER;k++) 
               {
                IN >> Value;
                if(Value!=0) *cR++ = Value;
                else 
                 {
                  VV = k;
                  for(k=k;k<MAXLAYER;k++) *cR++ = 0.0;
                  break;
                 }
               }
               if(k==MAXLAYER+1&&VV==0) 
                {
                    cerr << "# Error: Input::ReadParameters, " << Temp->Data.Name
                         << " token must end with 0.0" << endl;
                    cerr << "MAXLAYER = " << MAXLAYER << endl;                         
                    exit(NOT_END_WITH_ZERO_ERROR);
                };
              break;

       case INTMATRIX :
               cI = (int*) Temp->Data.Address;
               for(k=0;k<=N_ATOMS;k++)
               {
                IN >> Value;
                if(Value!=0) *cI++ = (int) Value;
                else break;
               }
               if(k==N_ATOMS+1)
                {
                    cerr << "# Error: Input::ReadParameters, " << Temp->Data.Name
                         << " token must end with 0.0" << endl;
                    cerr << "N_ATOMS = " << N_ATOMS << endl;
                    exit(NOT_END_WITH_ZERO_ERROR);
                };
              break;
       case BOOLEAN : 
              int *iI= (int*) Temp->Data.Address;
              IN >> CAR;     
              if     (strcmp(CAR,"On" )==0)  *iI= 1;
              else if(strcmp(CAR,"Off")==0)  *iI= 0;
              else if(strcmp(CAR,"Yes")==0)  *iI= 1;
              else if(strcmp(CAR,"No" )==0)  *iI= 0;
              else if(strcmp(CAR,"1"  )==0)  *iI= 1;
              else if(strcmp(CAR,"0"  )==0)  *iI= 0;
              else 
               {
                cerr << "# Boolean error ("<< Temp->Data.Name <<") in: "
                     << FileName << endl;
                exit(BOOLEAN_ERROR);
               };
              break;          

      };
      NVars ++;
    };
    Temp = PointerToNextData();
       };
    };

    // Close the stream
    IN.close();
    if( !TIP )
    {
     // Close the last layer
     Temp = PointerToFirstData();
     for( i=0; i<NumberOfElements(); i++)
      {
        if( Temp->Data.Type != LAYERDEFINITION ) 
                            Temp = PointerToNextData();
      };

     Y = (LayerDefinition*) Temp->Data.AuxiliaryAddress;
     cY = (List<LayerDefinition>*) Temp->Data.Address;
     (cY)->Add( (LayerDefinition) *Y);

     // Close the last simulation
     Temp = PointerToFirstData();
     for( i=0; i<NumberOfElements(); i++)
      {
        if( Temp->Data.Type != SIMULATIONDEFINITION ) 
                            Temp = PointerToNextData();
      };

     S = (SimulationDefinition*) Temp->Data.AuxiliaryAddress;
     cS = (List<SimulationDefinition>*) Temp->Data.Address;
     (cS)->Add( (SimulationDefinition) *S);
    };
    return NVars;
}

// Write parameters used to a file
int Input::WriteParameters()
{
    // Local variables
    int      NVars = 1;     // Number of variables writed
    LinkObject<Variable> *Temp; // Temporary variable

    // Open the stream
    ofstream OUT( FileName.c_str() );
    if (!OUT )
     {
      cerr << "Error: Input::WriteParameters, Can't open " << FileName << endl;
      return(-1);
     };

    // Write the parameters
    Temp = PointerToFirstData();
    while (NVars<=NumberOfElements())
    {
      if(Temp->Data.Type != XTALDEFINITION2 )
         OUT << Temp->Data;

      NVars ++;
      Temp = PointerToNextData();
    };

    // Close the stream
    OUT.close();
    return NVars;
}

// Output stream overloaded
ostream &operator<<( ostream &Out, Input A)
 {
  Out << " Filename : " << A.FileName << endl
      << (List<Variable>) A << endl;
  return Out;
 }
}
///////////////////////////////////////////////////////////////////////////
