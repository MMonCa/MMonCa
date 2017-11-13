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
/*
 *  body.h
 *  feliks
 *
 *  Created by Ignacio Romero on 11/9/06.
 *  Modified by I.R. on 31/5/07
 *
 * The body is one of the key ingredients of FELIKS. It is responsible for holding the nodes/elments, 
 * and each one of these  must belong to one and only one body. The power of such concept is that 
 * bodies can be activated or deactivated, duplicated,
 * they can interact pairwise, etc. It provides a layer of structure for elements and 
 * nodes that simplifies the treatment
 * of multibody problems, especially those with many interacting bodies with contact.
 *
 * The concept of body, as it is employed in this program, is not identical to that of "physical body". In fact,
 * if no bodies are defined in the input file, all nodes and elements will belong to a single "Default" body which
 * is always created, even if there are serveral "physical bodies" separated and without interaction. Similarly,
 * a body can be made of different materials, and actually, a body does not have a material associated with it.
 *
 * Also, body is a parent class for "regularbody". These are simple kind of bodies, with simple geometries, and made of
 * a single material. 
 */

#ifndef _body_h
#define _body_h


#include "Elements/element.h"
#include "Model/Node/node.h"
#include "Model/Parts/modelpart.h"
#include "Math/tensor.h"
#include "Model/Sets/faceset.h"


#include <map>
#include <vector>
#include <set>
#include <string>


namespace feliks
{
    namespace geometry
    {
        class BRep;
    }

    namespace meshing
    {
        class solidmesh;
    }
    
    namespace topology
    {
        class cell;
    }
}


class assembler;
class commandLine;
class energies;
class model;
class integrator;
class resultdata;


class body : public modelpart
{
	
public:
		
                                body(const std::string& name_, model& m, const size_t dim);
	virtual                     ~body();
	
	// scan the definition of the body from a commandline of the type 
	static body&                scan(const commandLine &cl, model &m);
	
	
	// initialization tasks
	virtual bool                check();
	virtual void                initialize();
	virtual void                parse();
	void                        parseInformationReset();
	void                        replaceNode(node& deleted, node& survivor);
	virtual void                setName(const std::string& newname);

	    
	virtual void                computeExplicitAcceleration();
    virtual void                computeExplicitEulerianAcceleration();
	virtual void                computeExplicitVelocity();
	virtual void                computeMass();
	const blue::ivector         computeMaxSurfaceVelocity() const;
	void                        computeReactionsSum(double *reac) const;
    virtual double              maxEigenvalueEstimate() const;

	
	// functions to retrieve components of the body
	int                         getCoordinationNumber(){return  coordinationNumber;}
    virtual feliks::meshing::solidmesh&         getTheSolidmesh();
    virtual const feliks::meshing::solidmesh&   getTheSolidmesh() const;
    virtual std::set<node*>     getNodesFrom0Cells(const std::set<feliks::topology::cell*>& cells) const;

	
	
	// functions to obtain information about the body
    virtual size_t              dimension() const;
	const blue::ivector         getCurrentCenterOfMass();
	double                      getMass() const;
	bool                        hasNode(const int label) const;

	
	// modify or update degrees of freedom of the body
	virtual void                advanceInTime(const integrator& i, const double dt);
	void                        incrementNCoord(){coordinationNumber++;}
	virtual void                incrementSolution(const integrator& i);
   	virtual void                remesh(bool& rec, bool& ghostnodes);
	void                        resetSolutionIncrements();
    void                        resetNodalMasses();
	void                        resetNodalForces();
	void                        rewindSolution();
	void                        setNCoord(const int c){coordinationNumber = c;}
	void                        setNodalMasses(const longvector& v);
    void                        updateCurrentState(const dofset::evaluation_time when);

    
    
	// display information or data
	virtual void                print(std::ostream &os=std::cout);

    std::map<std::string,faceset*>          theFacesets;


protected:

    size_t                      _dimension;
	bool                        hasVariableMass;
	int                         coordinationNumber;
    feliks::meshing::solidmesh*   theSolidmesh;

	double                      mass;
	bool                        massComputed;
    blue::istensor              inertia;
	
	void                        sortNodes();

    
private:

	// for the moment, we allow complete access to these classes
	friend class	assembler;
	friend class    gid;
	friend class    paraview;
	friend class    realtime;
	friend class    gmsh;
	friend class	asynchronousAnalysis;
	friend class	SPRsmoother;
	
	friend class	avi;
	friend class	backwardeuler;
	friend class	forwardeuler;
	friend class    hht;
	friend class	midpoint;
	friend class    newmark;
	friend class	centralDiff;
	friend class	integrator;
	friend class	hyperpara;
	friend class	quasiStatic;
	friend class	SCSCMatrix;
	friend class	CSCMatrix;
	friend class	petscmatrix;
    friend class    meshlessbody;
};

inline size_t      body::dimension() const {return _dimension;}


void	checkCoordNum(const double radio, std::vector<body*> &bodies);
void	printBodyElements(const body& b, std::ostream &of=std::cout);
void	printBodyFacesets(const body& b, std::ostream &of=std::cout);
void	printBodyNodeConnectivity(const body& b, std::ostream &of=std::cout);
void	printBodyNodes(const body& b, std::ostream &of=std::cout);
void	printBodyReactionsInConstraints(const body& b);


    
    
#endif

