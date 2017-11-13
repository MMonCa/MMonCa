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
#ifndef _INCLUDE_ATOMDEF_H_
#define _INCLUDE_ATOMDEF_H_
// Define Atom definition class ///////////////////////////////////////////

namespace BCA {

class AtomDefinition
{
 public:
    std::string Type;  // The chemical symbol
    int      Z;            // Atomic Number
    double     W;            // Atomic Weigth
    double     TAmp;         // Thermal RMS Displacement Amplitude
    double     TDebye;       // Debye temperature

    // Constructor
    AtomDefinition(const std::string &TC="H", double WC=1.0, int ZC=1, double TD=0.0 );
    // Compare two atoms
    int operator==( AtomDefinition &B );
    // Change the chemical symbol
    void SetName( const std::string &TC);
    // Get the chemical symbol
    const char* GetName();

    // Output streams overloaded : prototype
    friend std::ostream &operator<<(std::ostream &Out, AtomDefinition A);
    friend std::istream &operator>>(std::istream &In,  AtomDefinition &V );
};
}
#endif //_INCLUDE_ATOMDEF_H_
