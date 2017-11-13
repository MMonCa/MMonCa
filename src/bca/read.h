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
// READ.H
//
///////////////////////////////////////////////////////////////////////////
#ifndef _INCLUDE_READ_H_
#define _INCLUDE_READ_H_

// Includes
#include "layer.h"
#include "list.h"
#include "atom.h"
#include "crystal.h"
#include "simul.h"

#include <fstream>
#ifndef _toupper
#define _toupper(c) ((c) + 'A' - 'a')
#endif

// Defines

#define INTEGER               0
#define REAL                  1
#define CHAR                  2
#define VECTOR                3
#define LONG                  4
#define ATOMDEFINITION        5
#define XTALDEFINITION        6
#define LAYERDEFINITION       7
#define TEXT_                 8
#define ARRAY                 9
#define SIMULATIONDEFINITION 10
#define MATRIX               11
#define XTALDEFINITION2      12
#define INTMATRIX            13
#define BOOLEAN              14
#define MATRIX2              15

///////////////////////////////////////////////////////////////////////////

namespace BCA {

class Variable
{
   public:
    std::string   Name;
    char   Type;
    void   *Address;        // Variable address
    void   *AuxiliaryAddress; 

    // CONSTRUCTOR
    Variable( const std::string &N="", int T=REAL, void* AD=NULL, void* AA=NULL ) : Name(N)
     {
       Address          = AD;
       AuxiliaryAddress = AA;
       Type             = (char) T;
     };
    
    // Output stream overloaded
    friend std::ostream &operator<<(std::ostream &Out, Variable A);
};

///////////////////////////////////////////////////////////////////////////

class Input : public List<Variable>
{
   public:
	std::string FileName;
    int  CaseSensitive;

     // Constructor
    Input(const std::string &Name = "input_file", int CS = 0 ) : FileName(Name)
     {
       CaseSensitive = CS;
     };

    void Set(const std::string &Name = "input_file", int CS = 0 )
     {
       FileName = Name;
       CaseSensitive = CS;
     };

    // Main Method
    int ReadParameters( int TIP = 0 );
    int WriteParameters();

    // Output stream overloaded
    friend std::ostream &operator<<(std::ostream &Out, Input A);
};
}
#endif  // _INCLUDE_READ_H_
