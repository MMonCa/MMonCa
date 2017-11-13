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
/* node.cpp
*
* ignacio romero
* may 2000, revised sep 2003
* converted to c++ in jan 2006
*
* each node has associated certain number of local dofs, numbered starting from 0
*/

#include <iostream>
#include <iomanip>
#include <ostream>
#include <stdexcept>

#include <cstdlib>
#include <cstdarg>
#include <cmath>
#include <pthread.h>

#include "assert.h"
#include "Io/logger.h"
#include "Model/Node/node.h"
#include "General/idata.h"
#include "Io/message.h"
#include "Main/feliks.h"

#include "Math/Topology/topology.h"
#include "Math/tensor.h"
#include "Math/vector.h"

#ifdef WITHTBB
#include "tbb/spin_mutex.h"
#include "tbb/atomic.h"
#include "/usr/llvm-gcc-4.2/lib/gcc/i686-apple-darwin11/4.2.1/include/omp.h"
#endif

using namespace std;
using namespace blue;

int    node ::  maxlabel = -1;		// largest node label that has been employed


node :: node() :
	refCoords(),
	label(0), 
	invmass(0.0), damping(0.0), invdamping(0.0),
	constrDofs(0), hasLoad_(false), isOnBoundary_(false),
	glueSlave_(false), 
	masterGlueNode_(0),
	penetrationCount_(0),
	ndofs(0),			// dofs
	nodeTime(0.0)						// nodal forces
#ifdef WITHTBB
	,assemble_mutex()
#endif
{
    mass = 0.0;
}


node :: node(unsigned int labelc) : 
	refCoords(),
	label(labelc),
	invmass(0.0), damping(0.0), invdamping(0.0),
	constrDofs(0),  hasLoad_(false), isOnBoundary_(false),
	glueSlave_(false), masterGlueNode_(0),
	penetrationCount_(0),
	ndofs(0),
	nodeTime(0.0)
#ifdef WITHTBB
    ,assemble_mutex()
#endif
{
    mass     = 0.0;
    maxlabel = (labelc > maxlabel) ? labelc : maxlabel;
}



node :: node(unsigned int labelc, const ivector& c, const feliks::topology::cell *the0Cell) : 
    theCell(the0Cell),
	label(labelc),
	refCoords(c),
	invmass(0.0), damping(0.0), invdamping(0.0),
	constrDofs(0),  hasLoad_(false), isOnBoundary_(false),
	glueSlave_(false), masterGlueNode_(0),
	penetrationCount_(0),
	ndofs(0),   		// dofs
	nodeTime(0.0)
#ifdef WITHTBB
,assemble_mutex()
#endif
{
    mass     = 0.0;
    maxlabel = (labelc > maxlabel) ? labelc : maxlabel;
}





// copy constructor. A new node is created with the same data but different label.
node :: node(const node &nd)
#ifdef WITHTBB
: assemble_mutex()
#endif
{
	label       = ++maxlabel;
	ndofs       = nd.ndofs;
	constrDofs  = nd.constrDofs;
	hasLoad_	= nd.hasLoad_;
	penetrationCount_ = nd.penetrationCount_;
	isOnBoundary_	= nd.isOnBoundary_;
	glueSlave_  = nd.glueSlave_;
	masterGlueNode_ = nd.masterGlueNode_;
	mass        = nd.mass;
	invmass		= nd.invmass;
	damping		= nd.damping;
	invdamping	= nd.invdamping;
	nodeTime    = nd.nodeTime;

	refCoords   = nd.refCoords;	
}





/* returns all the memory allocated in NewNode and in InitializeNode
*/
node :: ~node()
{
}


std::ostream&  operator<<(std::ostream &os, const node &nd)
{	
	const int ndm(3);

	os << "\n\n Node  " << nd.label << flush;
    if (nd.theCell == 0) 
        os << "\n Cell empty";
    else
    {
        os << "\n Cell ";
        nd.theCell->print(os);
    }
    
	os << "\n   Mass                      " << setw(9) << nd.mass;
	os << "\n   Damping                   " << setw(9) << nd.damping;
	os << "\n   On boundary               " << (nd.isOnBoundary() ? "Yes" : "No");
	os << "\n   ID map                    " ;  nd.printID(os);
	
	// coordinates
	os << "\n   Reference position   ->";
	for (int k=0; k< ndm; k++)    os << "   " << right << setw(9) <<
        scientific << setw((int)os.precision()+6) << nd.refCoords(k);

	os << "\n   DOFs                   " ; 	nd.printDOFs(os);

	return os;
}



ivector& node :: acceleration(const dofset::evaluation_time when)
{
	throw std::runtime_error("This node type has no acceleration");
}



void  node :: assembleForce(const ivector& f)
{
	throw std::runtime_error ("Node should not assemble vector forces");
}



void  node :: assembleForce(const double f)
{
	throw std::runtime_error ("Node should not assemble scalar forces");
}


// careful with this function. It is to be used by problems
// where the coordinates do not change in time (heat, electromagnetism, ...)
ivector node :: coordinates(const dofset::evaluation_time when)
{
	return refCoords;
}


const ivector node :: coordinates(const dofset::evaluation_time when) const
{
	return const_cast<node&>(*this).coordinates(when);
}


ivector& node :: displacement(const dofset::evaluation_time when)
{
	throw std::runtime_error("This node type has no displacement");
}



ivector& node :: force()
{
	throw std::runtime_error("This node type has no force");
}


// glues the node to node2, making dU actually point to 
// the displacements of nd2, and replicating the ID map
// of nd2
// 
void node :: glueTo(const node& nd2)
{
	/*
	delete [] dU;
	dU = nd2.getSolutionIncrements();
	ID = nd2.getID();
	
	hasLoad_        = nd2.hasLoad();
	glueSlave_      = true;
	masterGlueNode_ = const_cast<node*>(&nd2);
	 */
}



// increment the nodal damping after locking access to it
void node :: incrementDamping(const double inc)
{
#ifdef WITHTBB
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(assemble_mutex);
#endif
	damping += inc;
	if (damping > 0.0) invdamping = 1.0/damping;


#ifdef WITHTBB
	my_lock.release();
#endif
}


// increment the nodal mass after locking access to it
void node :: incrementMass(const double inc)
{
    
#ifdef WITHTBB
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(assemble_mutex);
#endif
    
	mass += inc;
	if (mass > 0.0) invmass = 1.0/mass;

#ifdef WITHTBB
	my_lock.release();
#endif
}




/* initializes the data corresponding to the degrees of freedom (just for static problems).
* The data for problems with rotations and dynamic problems is done elsewhere.
* the node data ndofs should have been defined before calling this function, probably
* from an element function.
*/
void node :: initializeDofs()
{    
    if (ndofs == 0)  return;
}



void node :: initializeRotationalDofs(const iquaternion& qrot)
{
	logger::mainlog << "\nInitialization of rotational dofs not implemented for node " << getLabel();
}



void node :: initializeRotationalDofs(const ivector& director)
{
	logger::mainlog << "\nInitialization of director dofs not implemented for node " << getLabel();
}




// this function is required for the sort in the following function
bool node :: nodepLess(const node *ndp1, const node *ndp2) 
{return *ndp1 < *ndp2;}



void  node :: printForce(ostream &of)
{
	of << "\n" << "Label(" << setw(4) << label << ") :";
}




// raten = 1 for velocities
//       = 2 for accelerations
void  node :: printRates(const int raten, ostream &of) const
{		
	// if echo to screen display message
	of << "\n" << setw(4) << label;

	if ( getNodetype() == _Unode)
	{
		for (int k=0; k<3; k++)
			of << scientific << setw((int) of.precision()+8) << dynamic_cast<const Unode&> (*this).getUDS().velocity(dofset::tn1)(k);
	}
}




/* dumps the reactions in the nodal degrees of freedom to a file.
	Caution: the nodal forces must have been assembled before.
	The stored nodal force is the result of Fext - Fint - Ma
*/
void node :: printReactions(FILE *fp)
{
	extern double global_macheps;
	
    for (int i=0; i<3; i++)
	{
		double reac = force()(i);
		reac = fabs(reac) > 1000.0*global_macheps ? reac : 0.0;

		fprintf(fp, "  % 11.5e", reac);
	}
}



void node :: setInitialRate(const int local, const double v)
{
	if ( getNodetype() == _Unode || getNodetype() == _URnode || getNodetype() == _UPnode)
	{
		Unode& un = dynamic_cast<Unode&>(*this);
        
        if ( !getUDS().isConstrained(local) )
        {
            un.getUDS().velocity(dofset::tn)(local) = un.getUDS().velocity(dofset::tna)(local) =
                un.getUDS().velocity(dofset::tn1)(local) = v;
        }
	}
}
	


/*  stores ndofs as the number of degrees of freedom of the node unless
the node already has this information. In this case only write ndofs
if it is larger that what already existed. This is done to find out
how many dofs a node has from the information coming from the connecting
elements. In general, a node must have as many dofs as the maximum
number of dofs needed by a connecting element
*/
void   node :: setNdofs(int n)
{
    ndofs =  n > ndofs ? n : ndofs;
}


void node :: setRotation(const irotation& r)
{
	//	if (rottype > ROT_NONE)	rotn1->extractFromRotation(r);
}



ivector& node ::	velocity(const dofset::evaluation_time when)
{
	throw std::runtime_error("This node type has no velocity");
}




Unode :: Unode(const int label_, const ivector& coor, const feliks::topology::cell *the0Cell) :
	node(label_, coor, the0Cell),
	_currCoordinates(coor),
	U()
{
}


 Unode :: Unode(const int label_) :
	node(label_),
	U()
{
}


void Unode ::	advanceIntegration(const integrator& integ, const double dt)
{
	integ.advanceDofset( U);
	_currCoordinates = refCoords + U.displacement();
}



void Unode :: assembleForce(const ivector& f)
{
#ifdef WITHTBB
	// each node locks its assemble
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(assemble_mutex);
	U.force() += f;
	my_lock.release();
    
#else
	U.force() += f;	
#endif
}



void Unode :: assembleForceDof(const int dof, const double f)
{
	assert(dof<3);
	U.force()(dof) += f;
}



// compute the acceleration, branching if there are constraints to make it faster
void Unode :: computeExplicitAcceleration()
{
	double im = getInvMass();
	if (im == 0.0)
		logger::mainlog << "\n Zero mass in node " << getLabel() << flush;
	
    // negative masses can happen in meshless solutions. Set ballistic motion
    extern double global_macheps;

    if (im <=0.0 || getMass() <= global_macheps)
    {
        cout << " mingetMass!\n ";
       U.acceleration(dofset::tna).setZero();
    }
    else    
    {
        if (constrDofs == 0)
        {
            U.acceleration() = U.force() * im;
        }
        
        else
        {
            // acceleration from forces
            if ( !U.isConstrained(0) ) U.acceleration(dofset::tna)(0) = U.force()(0) * im;
            if ( !U.isConstrained(1) ) U.acceleration(dofset::tna)(1) = U.force()(1) * im;
            if ( !U.isConstrained(2) ) U.acceleration(dofset::tna)(2) = U.force()(2) * im;
        }
    }
}


void Unode ::	computeExplicitEulerianAcceleration()
{	
}




void Unode ::	computeExplicitVelocity()
{	
	double idamp = getInvDamping();
	
	// velocity from forces
	if ( !U.isConstrained(0) ) U.velocity(dofset::tna)(0) = U.force()(0) * idamp;
	if ( !U.isConstrained(1) ) U.velocity(dofset::tna)(1) = U.force()(1) * idamp;
	if ( !U.isConstrained(2) ) U.velocity(dofset::tna)(2) = U.force()(2) * idamp;	
}



void Unode :: constrainDof(const int local)
{
	assert(local<3);
	if ( U.setConstrainedDof(local) ) constrDofs++;
}


ivector Unode :: coordinates(const dofset::evaluation_time when)
{
	if  (when == dofset::tna)
        return getReferenceCoordinates() + U.displacement(dofset::tna);
		
	else if (when == dofset::t0) 	
		return refCoords;
	
	else if	(when == dofset::tn)
			return getReferenceCoordinates() + U.displacement(dofset::tn);
		
	else if (when == dofset::tn1)
			return getReferenceCoordinates() + U.displacement(dofset::tn1);
		
	else 
		return _currCoordinates;
}


const ivector Unode :: coordinates(const dofset::evaluation_time when) const
{
	return const_cast<Unode&>(*this).coordinates(when);
}





void Unode :: incrementSolution(const integrator& integ)
{
	integ.incrementDofset(U);
	_currCoordinates = refCoords + U.displacement();
}




void Unode :: localizeSolutionIncrement(const longvector& delta,  partitionDofs pt)
{
	U.localizeSolutionIncrement(delta);
}



void  Unode :: printDOFs(ostream &of) const
{
	U.printDOFS(of);
}



void  Unode :: printForce(ostream &of)
{
	node::printForce(of);
	getUDS().force().print(of);
}

void Unode :: recoverFromBackupDofs()
{
	U.recoverFromLastSolution();
}



void Unode :: setDeltaDofsToZero()
{
	U.resetDelta();
}


void Unode :: setDof(const int local, const double v)
{
	assert(local<3);
	U.setDof(local, v);
}


void Unode :: setDofIncrement(const int local, const double v)
{
	assert(local<3);
	U.increment()(local) = v;
}

void Unode ::	setInitialDof(const int local, const double v)
{
	assert(local<3);
	U.displacement(dofset::tn)(local) = U.displacement(dofset::tn1)(local) = U.displacement()(local) = v;
}





UPnode :: UPnode(const int label_, const ivector& coor, const feliks::topology::cell *the0Cell) :
	Unode(label_, coor, the0Cell)
{
}



void UPnode ::	advanceIntegration(const integrator& integ, const double dt)
{
	Unode::advanceIntegration(integ, dt);
	integ.advanceDofset( P );
}



void UPnode :: assembleForce(const double f)
{
#ifdef WITHTBB
	// each node locks its assemble
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(assemble_mutex);
	
	P.force() += f;
	
	my_lock.release();
	
	/*
	 // one mutex lock per node. Many assemblies can be performed simultaneously
	 // what must be avoided is assembly over the same node
	 pthread_mutex_lock( &assemble_mutex);
	 for (int i=0; i< ndofs; i++) force[i] += f[i];
	 pthread_mutex_unlock( &assemble_mutex);
	 */
#else
	P.force() += f;
#endif
}


void UPnode :: assembleForceDof(const int dof, const double f)
{
	assert(dof<4);
	
	if (dof < 3) U.force()(dof) += f;
	else		 P.force()	  += f;
}



void UPnode :: constrainDof(const int local)
{
	if (local < 3)  Unode::constrainDof(local);
	else if(P.setConstrainedDof(local-3) ) constrDofs++;
}



void UPnode :: incrementSolution(const integrator& integ)
{
	Unode::incrementSolution(integ);
	integ.incrementDofset(P);
}




void UPnode :: localizeSolutionIncrement(const longvector& delta,  partitionDofs pt)
{
	if (pt == PARTITION_1 || pt == PARTITION_ALL)	U.localizeSolutionIncrement(delta);
	if (pt == PARTITION_2 || pt == PARTITION_ALL) 	P.localizeSolutionIncrement(delta);
}



void  UPnode :: printDOFs(ostream &of) const
{
	U.printDOFS(of);
	P.printDOFS(of);
}


void  UPnode :: printForce(ostream &of)
{
	node::printForce(of);
	getUDS().force().print(of);
	of << getPDS().force();
}


void UPnode :: recoverFromBackupDofs()
{
	U.recoverFromLastSolution();
	P.recoverFromLastSolution();
}



void UPnode :: setDeltaDofsToZero()
{
	Unode::setDeltaDofsToZero();
	P.resetDelta();
}


void UPnode ::	setDof(const int local, const double v)
{
	assert(local < 4);
	if (local < 3)  U.setDof(local, v);
	else			P.setDof(v);
}


void UPnode ::	setDofIncrement(const int local, const double v)
{
	assert(local < 4);
	if (local < 3)  U.increment()(local) = v;
	else			P.increment() = v;
}



void UPnode ::	setInitialDof(const int local, const double v)
{
	assert(local<4);
	if (local<3)
		U.displacement(dofset::tn)(local) = U.displacement(dofset::tn1)(local) = U.displacement()(local) = v;
	else
		P.displacement(dofset::tn) = P.displacement(dofset::tn1) = P.displacement() = v;
}




UDnode :: UDnode(const int label_, const ivector& coor, const feliks::topology::cell *the0Cell) :
	Unode(label_, coor, the0Cell)
{
}


void UDnode ::	advanceIntegration(const integrator& integ, const double dt)
{
	Unode::advanceIntegration(integ, dt);
	integ.advanceDofset( D );
}


void UDnode :: assembleForceDof(const int dof, const double f)
{
	assert(dof<5);
	
	if (dof < 3) U.force()(dof)   += f;
	else		 D.force()[dof-3] += f;
}



void UDnode :: constrainDof(const int local)
{
	if (local < 3)  Unode::constrainDof(local);
	else if(D.setConstrainedDof(local-3) ) constrDofs++;
}


void UDnode :: incrementSolution(const integrator& integ)
{
	Unode :: incrementSolution(integ);
	integ.incrementDofset(D);
}




// initialize rotational dofs with the rotation that takes E3 to director
// with no drill. If the rotational dofs are already initialized, return
void UDnode :: initializeRotationalDofs(const ivector& director)
{	
	// shell directors
	irotation rot;
	rot.beRotationWithoutDrill(director);
	
	getDDS().rotation(dofset::t0).extractFromRotation(rot);
	getDDS().rotation(dofset::tn) = getDDS().rotation(dofset::tna) = getDDS().rotation(dofset::tn1) =
        getDDS().rotation(dofset::t0);
}




void UDnode :: localizeSolutionIncrement(const longvector& delta,  partitionDofs pt)
{
	if (pt == PARTITION_1 || pt == PARTITION_ALL) 	U.localizeSolutionIncrement(delta);
	if (pt == PARTITION_2 || pt == PARTITION_ALL)	D.localizeSolutionIncrement(delta);
}



void  UDnode :: printDOFs(ostream &of) const
{
	U.printDOFS(of);
	D.printDOFS(of);
}


void  UDnode :: printForce(ostream &of)
{
	node::printForce(of);
	getUDS().force().print(of);
	getDDS().force().print(of);
}


void UDnode :: recoverFromBackupDofs()
{
	U.recoverFromLastSolution();
	D.recoverFromLastSolution();
}



void UDnode :: setDeltaDofsToZero()
{
	Unode::setDeltaDofsToZero();
	D.resetDelta();
}


void UDnode ::	setDofIncrement(const int local, const double v)
{
	assert(local < 5);
	if (local < 3)  U.increment()(local) = v;
	else			D.increment()[local-3] = v;
}



void UDnode ::	setInitialDof(const int local, const double v)
{
	assert(local<5);
	if ( local < 3) 
		U.displacement(dofset::tn)(local) = U.displacement(dofset::tn1)(local) = U.displacement()(local) = v;
	else
		;
}
		




Pnode :: Pnode(const int label_, const ivector& coor, const feliks::topology::cell *the0Cell) :
	node(label_, coor, the0Cell)
{
}




void Pnode ::	advanceIntegration(const integrator& integ, const double dt)
{
	integ.advanceDofset( P );
}



void Pnode :: assembleForce(const double f)
{
#ifdef WITHTBB
	// each node locks its assemble
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(assemble_mutex);
	
	P.force() += f;
	
	my_lock.release();
	
	/*
	 // one mutex lock per node. Many assemblies can be performed simultaneously
	 // what must be avoided is assembly over the same node
	 pthread_mutex_lock( &assemble_mutex);
	 for (int i=0; i< ndofs; i++) force[i] += f[i];
	 pthread_mutex_unlock( &assemble_mutex);
	 */
#else
	P.force() += f;
#endif
}


void Pnode :: assembleForceDof(const int dof, const double f)
{
	assert(dof<1);
	P.force() += f;
}



void Pnode ::	computeExplicitAcceleration()
{
	double im = getInvMass();
	
	// acceleration from forces
	if ( !P.isConstrained(0) ) P.acceleration(dofset::tna)= P.force() * im;
}


void Pnode ::	computeExplicitVelocity()
{	
	double idamp = getInvDamping();
	
	// velocity from forces
	if ( !P.isConstrained(0) ) P.velocity(dofset::tna) = P.force() * idamp;
}


void Pnode :: constrainDof(const int local)
{
	if(P.setConstrainedDof(0) ) constrDofs++;
}



void Pnode :: incrementSolution(const integrator& integ)
{
	integ.incrementDofset(P);
}




void Pnode :: localizeSolutionIncrement(const longvector& delta,  partitionDofs pt)
{
	P.localizeSolutionIncrement(delta);
}




void  Pnode :: printDOFs(ostream &of) const
{
	P.printDOFS(of);
}


void  Pnode :: printForce(ostream &of)
{
	node::printForce(of);
	of << getPDS().force();
}


void Pnode :: recoverFromBackupDofs()
{
	P.recoverFromLastSolution();
}



void Pnode :: setDeltaDofsToZero()
{
	P.resetDelta();
}




void Pnode :: setDof(const double v)
{
	P.setDof(v);
}


void Pnode ::	setDofIncrement(const int local, const double v)
{
	P.increment() = v;
}


void Pnode ::	setInitialDof(const int local, const double v)
{
	assert(local<1);
	P.displacement(dofset::tn) = P.displacement(dofset::tn1) = P.displacement() = v;
}





template <int N>
Vnode<N> :: Vnode(const int label_, const ivector& coor, const feliks::topology::cell *the0Cell) :
node(label_, coor, the0Cell)
{
}




template <int N>
void Vnode<N> :: advanceIntegration(const integrator& integ, const double dt)
{
	integ.advanceDofset( ds, dt);
}



template <int N>
void Vnode<N> :: assembleForce(const double* f)
{
#ifdef WITHTBB
	// each node locks its assemble
	tbb::spin_mutex::scoped_lock my_lock;
	my_lock.acquire(assemble_mutex);
	
	for (int i=0; i<N; i++) ds.force()[i] += f[i];
	
	my_lock.release();
	
	/*
	 // one mutex lock per node. Many assemblies can be performed simultaneously
	 // what must be avoided is assembly over the same node
	 pthread_mutex_lock( &assemble_mutex);
	 for (int i=0; i< ndofs; i++) force[i] += f[i];
	 pthread_mutex_unlock( &assemble_mutex);
	 */
#else
	for (int i=0; i<N; i++) ds.force()[i] += f[i];
#endif
}


template <int N>
void Vnode<N> :: assembleForceDof(const int dof, const double f)
{
	assert(dof<N);
	ds.force()[dof] += f;
}


template <int N>
void Vnode<N> ::	computeExplicitAcceleration()
{
	double im = getInvMass();
	
	// acceleration from forces
	for (int i=0; i<N; i++)
		if ( !ds.isConstrained(i) ) ds.acceleration(dofset::tna)[i] = ds.force()[i] * im;
}


template <int N>
void Vnode<N> ::computeExplicitVelocity()
{	
	double idamp = getInvDamping();
	
	// velocity from forces
	for (int k=0; k<N; k++)
		if ( !ds.isConstrained(k) ) ds.velocity(dofset::tna)[k] = ds.force()[k] * idamp;
}


template <int N>
void Vnode<N> :: constrainDof(const int local)
{
	if(ds.setConstrainedDof(local) ) constrDofs++;
}


template <int N>
void Vnode<N> :: incrementSolution(const integrator& integ)
{
	integ.incrementDofset(ds);
}



template <int N>
void Vnode<N> :: localizeSolutionIncrement(const longvector& delta,  partitionDofs pt)
{
	ds.localizeSolutionIncrement(delta);
}



template <int N>
void  Vnode<N> :: printDOFs(ostream &of) const
{
	ds.printDOFS(of);
}


template <int N>
void  Vnode<N> :: printForce(ostream &of)
{
	node::printForce(of);
	getVDS().printForce(of);
}

template <int N>
void Vnode<N> :: recoverFromBackupDofs()
{
	ds.recoverFromLastSolution();
}


template <int N>
void Vnode<N> :: setDeltaDofsToZero()
{
	ds.resetDelta();
}


template <int N>
void Vnode<N> :: setDofIncrement(const int local, const double v)
{
	assert(local < N);
	ds.increment()[local] = v;
}


template <int N>
void Vnode<N> :: setInitialDof(const int local, const double v)
{
	assert(local<N);
	ds.displacement(dofset::tn)[local] = ds.displacement(dofset::tn1)[local] = ds.displacement()[local] = v;
}





std::ostream&  operator<<(std::ostream &os, const Unode &nd)
{	
	os << dynamic_cast<const node&>(nd) << "\n";
	os << nd.getUDS() << flush;
	
	return os;
}

/*

std::ostream&  operator<<(std::ostream &os, const URnode &nd)
{	
	os << dynamic_cast<const node&>(nd) << endl;
	os << nd.getUDS() << endl;
	os << nd.getRDS() << flush;
	
	return os;
}
*/


std::ostream&  operator<<(std::ostream &os, const UDnode &nd)
{	
	os << dynamic_cast<const node&>(nd) << "\n";
	os << nd.getUDS() << "\n";
	os << nd.getDDS();
	
	return os;
}


std::ostream&  operator<<(std::ostream &os, const Pnode &nd)
{	
	os << dynamic_cast<const node&>(nd) << "\n";
	os << nd.getPDS();
	
	return os;
}


template <int N>
std::ostream&  operator<<(std::ostream &os, const Vnode<N> &nd)
{	
	os << dynamic_cast<const node&>(nd) << "\n";
	os << nd.getVDS();
	
	return os;
}



