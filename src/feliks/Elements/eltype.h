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
/* eltype.h
 *
 * iro, november 2006
 *
    As provided by the user input, a table of all the elements types is created and stored in
	memory, assigning a label to each type. This table is called allEltypes and is created
	once for all at the very beginning of the program execution (and deleted at the very end).
 
 * in this function the element types are defined. Objects in this class serve two purposes.
 *
 * 1) 
 */

#pragma once
#ifndef _feliks_eltype_h
#define _feliks_eltype_h

#include <vector>
#include <string>
#include <map>
#include <iostream>

#include "Analysis/energies.h"
#include "Math/tensor.h"
#include "Model/Node/node.h"
#include "Io/usercommand.h"


namespace feliks
{
    class elementState;
    namespace topology
    {
        class cell;
    }
}



class element;
class evalspot;
class body;
class shapefun;

/* type of element geometries. This is used for plotting purposes. For example: a 3 node element
   can be a triangle, a bar, or a spring. Also, a 4 node element can be a bar, a tetrahedra,
   a shell, a contact element ... There is no way of distinguishing among
   them from the node list only, thus a "geometry" must be defined too. This is implicitely done when
   chosing the eltype, hence the "geometry" is stored in the eltype structure */
enum geometryT{
	NOPLOT =-1 ,
	POINT  ,
	CURVE  ,
	SURFACE,
	VOLUME,
    COHESIVE
};




class eltype
{
public:

                            eltype(commandLine &cl);
                            eltype();
    virtual                 ~eltype();
	void                    addToList();
    static std::map<int, eltype*> allEltypes;

	
    
	// Operations on the global list of element types
	static void             freeEltypes();
	static eltype&          getEltypeFromLabel(int label);
	static size_t              getNEltypes();
	static void             listEltypes(ostream &of=std::cout);
	static void             scan(const commandLine &cl);
	
	
	// functions to get some data from one eltype
	int                     getDofsPerNode() const;
	const geometryT&		geometry() const;
	const nodetypeT&        nodetype() const;

	const int&              label()	const;
    const int&              materialLabel() const;
    const std::string&           name() const;

	// information
	bool                    isVisible() const;
	virtual void            print(ostream &of=cout);

	
	
	// this is the main creation function, which must be specified in each element type
	virtual element*        createElement(const int lbl, const feliks::topology::cell*, const int nv, node **ndlist, feliks::elementState& es) const=0;
    virtual evalspot*       createEvalspot(const double volume, const std::vector<shapefun>& theShp, const double initialspacing) const;
    virtual element*        createRandomElement() const {return 0;} // XX needs to be done

	
	// this function creates a node of the type associated with the eltype
	node*                   createNode(const int label, const blue::ivector& coor) const ;
    node*                   createNode(const int label, const blue::ivector& coor, const feliks::topology::cell* the0Cell) const;
    
	

    std::vector<double>     geoDimensions;				// geometric dimensions and other possible data
	int                     _label;

protected:
	
	// a table with all the user-defined element types
	static int              largestEltypeLabel;
	
	// description of the element type and basic data
	std::string                 _name;
	std::string                 formulation;					// formulation type
	std::string                 explanation;					// a brief explanation of the element
	int                         _materialLabel;					// label of the material of the elementtype
	geometryT                   _geometry;                      // for plotting. See above
	nodetypeT                   _nodetype;
	
    std::vector<std::string>    geoNames;					// the names of these data
	
	friend class element;
};

inline const int&           eltype :: label()	const               {return _label;}
inline const geometryT&     eltype :: geometry() const              {return _geometry;}
inline const int&           eltype :: materialLabel() const         {return _materialLabel;}
inline const std::string&   eltype :: name() const                  {return _name;}
inline const nodetypeT&     eltype :: nodetype() const              {return _nodetype;}
inline bool                 eltype :: isVisible() const             {return _geometry != NOPLOT;}





#endif
