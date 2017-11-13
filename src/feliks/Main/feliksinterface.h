/* Finite Element Method Module
 *
 * Author: ignacio.romero@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain,
 *      and       Technical University of Madrid (UPM), Madrid, Spain
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
 */ 
//
//  feliksinterface.h
//
//
//  Created by Ignacio Romero on 3/5/13.
//  Modified by I. Martin-Bragado. Feb 2015
//


#ifndef _feliksinterface_h
#define _feliksinterface_h

#include <vector>
#include <map>
#include <string>
#include "io/Diagnostic.h"

namespace feliks
{
    
    
    // material description
    // this description defines the mechanical response of a pure phase
    // the behavior of alloys is given, at this moment, by a weighted sum of responses
    // of the phases
    class mmoncaMaterial
    {
    public:
        mmoncaMaterial(const std::string& name,
                       double thermalexp,
                       double ref_temperature,
                       double E,
                       double nu,
                                  //hash            alpha    beta
                       std::map<std::string, std::pair<double, double> > eigenstrains
                       ) :
        name_(name), thexpansion_(thermalexp), ref_temp_(ref_temperature), E_(E), nu_(nu),
        	eigenstrains_(eigenstrains) {}
        
        std::string             name_;
        double                  thexpansion_;            // thermal expansion coeff
        double                  ref_temp_;               // reference temperature
        double                  E_, nu_;                 // moduli
        std::map<std::string, std::pair<double, double> > eigenstrains_;
    };
    
    
    
    
    // the state of a finite element
    // at the very instant of creation, the number of possible phases must be
    // defined and never modified. If desiered, the composition can then be
    // set to 1 for one phase, and 0 for the rest, but never altered
    class elementState
    {
    public:
        elementState(const int nmaterials)
        {
            composition_.resize(nmaterials);
        }
        
        
        
        // this is the state of the element.
        // it must be set before the FE computation. The composition (in unary fractions) may
        // change from step to step. The temperature and
        // initial strains are identical for all phases
        void setState(const std::vector<double>& composition,
                      double temp,
                      double initx,
                      double inity,
                      double initz,
                      const std::map<std::string, unsigned> &defects)
        {
            if ( composition.size() != composition_.size() )
                WARNINGMSG("Warning in elementState assignment");
            
            composition_ = composition;
            temperature_ = temp;
            iexx_        = initx;
            ieyy_        = inity;
            iezz_        = initz;
            defects_     = defects;
        }
        
        
        std::vector< double >   composition_;                   // fraction of each of the materials
        double                  temperature_;
        double                  iexx_, ieyy_, iezz_;            // initial strains
        double                  sxx, syy, szz, syz, sxz, sxy;   // stress
        double                  exx, eyy, ezz, eyz, exz, exy;   // strain
        std::map<std::string, unsigned> defects_;	            // number of defects
    };
    
    class mmoncaInterface
    {
        
    public:
        // the constructor of the interface allows to share information that is never
        // going to change (vertices, cells, boundc, materials) or whose memory position
        // is never going to change (element state). This makes transfer of data
        // between codes much easier and faster
        mmoncaInterface(const std::vector<double>&          vertices,
                        const std::vector<int>&             cells,
                        std::vector<elementState>&          es,
                        const std::map<int,char>&           boundc,
                        const std::vector<mmoncaMaterial>&  mat)
        :
        vcoordinates_(vertices),
        connectivity_(cells),
        bc_(boundc),
        theMaterials_(mat),
        theElements_(es)
        {
            disp_.resize(vcoordinates_.size());
        }
        
        // return the displacements after computation
        std::vector<double>&  displacements() { return disp_;}
        
        
        // vector of dimension 3*N, N: number of vertices.
        // it holds the coordinates of the vertices
        // in the mesh in a consecutive fashion: X1, Y1, Z1, X2, Y2, Z2, ...
        const std::vector< double >&      vcoordinates_;
        
        // vector of dimension 8*M, M: number o hexahedra
        // it holds the labels (pointers) of the vertices of each hexahedron
        // note: first vertex has label 0
        const std::vector< int >&         connectivity_;
        
        // map holding the boundary conditions of the vertices which are not
        // free (if a vertex has no boundary conditions, then don't insert
        // anything in the map
        // bc(i,j) holds the boundary condition for the i-vertex
        // the mask j holds all the conditions:
        //      j & 1, if x dof is constrained
        //      j & 2, if y dof is constrained
        //      j & 4, if z dof is constrained
        const std::map< int, char >&      bc_;
        
        
        // the array of materials
        //
        const std::vector<mmoncaMaterial>&      theMaterials_;
        
        
        // all the information on the elements. This might be changed
        // to re-use the coordinates, connectivity and bc
        std::vector< elementState >&     theElements_;
        
        
        // vector holding the vertex displacements after the computations
        // the dimension is 3*N, and the data is stored in consecutive
        // fashion as in U1, V1, W1, U2, V2, W2, U3 ...
        std::vector< double >           disp_;
        
        ;
    };
}

#endif
