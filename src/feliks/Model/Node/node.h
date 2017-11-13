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
/* node.h
 *
 * ignacio romero
 * date: 4 25 00, revised august 2003
 * converted to C++ class in january 2006
 * Unode, UDnode, URnode subclasses created in july 2009
 *
 * the node is the data structure that holds the information of the geometric definition of
 * a model vertex as well as all the solution data in a FE node (dofs, rates, etc).
 *
 *
 * Usage: the first thing one must do, even before defining any node is setting the
 * spatial dimension of the model they belong to. This is done wih the funcion
 * setNodesDimension. After this, nodes can be created and deleted at will.
 *
 * After all the nodes have been created, elements should call their corresponding
 * nodes telling them the number of dofs with the funcion setNDofs(). In this way,
 * the program takes the ndofs in one node as the maximum number of dofs/node in
 * all connecting elements.
 *
 * Next, space for the dofs is created with the function initializeDofs(), but this
 * must be done in the sequence just explained, because the memory allocator needs
 * information about the number of dofs per node.
 *
 */

#ifndef _node_h
#define _node_h

	
#include <iostream>
#include <ostream>
#include <vector>
#include "Io/Postprocessor/resultdata.h"
#include "Math/vector.h"
#include "Math/tensor.h"
#include "Analysis/Dofs/dofset.h"

#ifdef WITHTBB
#include "tbb/spin_mutex.h"
#endif

namespace feliks
{
    namespace topology
    {
        class cell;
    }
}


class element;
class integrator;



using namespace std;

	


enum	nodetypeT{_EMPTYnode, _Pnode, _Unode, _UPnode, _URnode, _UDnode, _SPnode, _Vnode};


class node
{
	
public:
	
	int                         label;			// node label. Not necessarily consecutive (must be unique)
    size_t                      ndofs;			// total number of dofs in node (active and constrained)
	resultdata                  theSmoothedResult;
		
                                node();
                                node(unsigned int label_);
                                node(unsigned int label_, const blue::ivector& c, const feliks::topology::cell *the0Cell=0);
                                node(const node &nd);
	virtual                     ~node();
	
		
	// integration tasks
	virtual	void                advanceIntegration(const integrator& integ, const double dt)=0;
	virtual	void                incrementSolution(const integrator& integ)=0;
	virtual void                computeExplicitAcceleration()=0;
	virtual void                computeExplicitVelocity()=0;
	        virtual void            computeExplicitEulerianAcceleration()=0; 
	void                        glueTo(const node& nd2);
	void                        initializeDofs();
	virtual void                initializeRotationalDofs(const blue::iquaternion& q);
	virtual void                initializeRotationalDofs(const blue::ivector& director);
	inline void                 resetPenetrationCount(){penetrationCount_=0;}
	inline void                 setHasLoad()					{hasLoad_ = true;}
	inline void                 setOnBoundary()					{isOnBoundary_=true;}
	void                        setNdofs(const int n);
	void                        setRotation(const blue::irotation& r);
	inline void                 setLabel(const int k){label = k;}
	inline void                 setTime(const double t)	{nodeTime = t;}
	
	
	// Connectivity data. Elements touching the node
	void                        setConnectivity(const std::vector<element*> inElements_);

	
	// functions to get information about the node
    inline const feliks::topology::cell*    getCellP() const     {return theCell;}
	inline int                  getLabel() const				{return label;}
	virtual	int                 getNDofs() const=0;
	int                         getNOrphanDofs() const;
	inline node&                getGlueMaster()					{return *masterGlueNode_;}
	virtual nodetypeT           getNodetype() const=0;
	inline int                  getPenetrationCount() const		{return penetrationCount_;}
	inline double               getTime() const					{return nodeTime;}
	inline bool                 hasLoad() const					{return hasLoad_;}
	inline bool                 isGlued() const					{return glueSlave_;}
	inline bool                 isOnBoundary() const			{return isOnBoundary_;}
	virtual bool                isOrphan() const=0;
	
	// get degrees of freedom and time derivatives
	virtual			blue::ivector&	acceleration(const dofset::evaluation_time when=dofset::tna);
	virtual			blue::ivector		coordinates(const dofset::evaluation_time when=dofset::tna);
	virtual const	blue::ivector		coordinates(const dofset::evaluation_time when=dofset::tna) const;
	virtual			blue::ivector&	displacement(const dofset::evaluation_time when=dofset::tna);
	inline const	blue::ivector&	getReferenceCoordinates() const	{return refCoords;}
	virtual			blue::ivector&	velocity(const dofset::evaluation_time when=dofset::tna);
	
	// functions to get dofsets. Declare all of them, irrespective of the
	// node type so that dynamic_casting is avoided in some function calls
	virtual vectorDofset&		getUDS()    {cout << "should not be called, getUDS"; return *(new vectorDofset);}
	virtual const vectorDofset& getUDS() const {cout << "should not be called, getUDS"; return *(new vectorDofset);}
	
	
	// functions to manipulate dofs
	virtual void                recoverFromBackupDofs()=0;
	inline void                 incrementPenetrationCount(){penetrationCount_++;}
	void                        setDofs(longvector dofs);
	virtual void                setDofIncrement(const int local, const double v)=0;
	virtual void                setInitialDof(const int local, const double v)=0;
	void                        setInitialRate(const int local, const double v);
	virtual void                setDeltaDofsToZero()=0;
	virtual void                setReferenceCoordinates(const blue::ivector& c){refCoords = c;}

	
	// functions to assign and retrieve dof/equation mapping ID. Local id start from 0
	virtual void                constrainDof(const int local)=0;
	virtual int                 getID(const int local) const=0;

	
	// functions to retrieve information about free/constrained dofs
	virtual int                 getNConstrainedDofs() const=0;
	inline short                getNFreeDofs() const					{return ndofs - constrDofs;}
	inline bool                 hasConstrainedDofs() const				{return (constrDofs > 0);}
	virtual bool                isDofConstrained(const int local) const=0;
	

	// functions to manipulate the nodal mass
	void				dampingReset();
	void				resetMass();
    void                changeMass(double m);
    
	double				getDamping() const;
	double				getMass() const;
	double				getInvDamping() const;
	double				getInvMass() const;
	void				incrementDamping(const double inc);
	void				incrementMass(const double inc);
	
	
	// functions to manipulate the forces
	virtual void		assembleForceDof(const int dof, const double f)=0;
	virtual	void		assembleForce(const blue::ivector& f);
	virtual	void		assembleForce(const double   f);
	virtual blue::ivector&	force();
	virtual void		forceReset()=0;

	virtual void		localizeSolutionIncrement(const longvector& delta,  partitionDofs pt=PARTITION_ALL)=0;


	// functions to display or print information about the node
	// and its dofs (1,2,3,...)
	virtual void		printID(ostream &of=cout) const=0;
	virtual void		printDOFs(ostream &of=cout) const =0;
	virtual void		printForce(ostream &of=cout);
	void				printRates(const int raten, ostream &of=cout) const;
	void				printReactions(FILE *fp=stdout);


	// this function defines an order relation by comparing node labels
	friend std::ostream&  operator<<(std::ostream &os, const node &nd);
	inline friend bool operator<(const node &nd1, const node &nd2) {return nd1.label < nd2.label;}
	inline friend bool operator==(const node &nd1, const node &nd2) {return nd1.label == nd2.label;}
	static bool nodepLess(const node *ndp1, const node *ndp2);
	
	
private:
	static int              maxlabel;		// largest node label that has been employed (may be <0)
	
    const feliks::topology::cell    *theCell;

    double                  mass;
	double                  invmass;        // lumped nodal mass
	double                  damping;
    double                  invdamping;
	double                  nodeTime;		// for integrators with individual times
	bool                    hasLoad_;
	bool                    isOnBoundary_;	// true if the node is on the boundary of the model
	bool                    glueSlave_;		// true if glued to another node from which it obtains dofs
	int                     penetrationCount_;
	node*                   masterGlueNode_;
	
	node& operator=(const node &nd);	// avoid copying nodes
	
	
protected:
	size_t				constrDofs;		// number of constrained dofs
	blue::ivector				refCoords;		// reference coordinates of the node
    
#ifdef WITHTBB
	tbb::spin_mutex		assemble_mutex;
#endif

#ifdef WITHMPI
public:
	vector<int>			neighborPartitions;
#endif
	
	
};

inline void		node :: dampingReset()			{damping = invdamping = 0.0;}
inline void		node :: resetMass()				{mass = invmass = 0.0;}

inline double	node ::	getInvDamping() const	{return invdamping;}
inline double	node ::	getInvMass() const      {return invmass;}
inline double	node :: getDamping() const		{return damping;}
inline double	node :: getMass() const			{return mass;}
inline void     node :: changeMass(double m) 	{mass = m;}







class Unode : public node
	{
		
	public:
		Unode(const int label_, const blue::ivector& coor, const feliks::topology::cell *the0Cell=0);
		Unode(const int label_);
		virtual ~Unode(){}
		
		
		blue::ivector&                acceleration(const dofset::evaluation_time when=dofset::tna){return U.acceleration(when);}
		virtual	void            advanceIntegration(const integrator& integ, const double dt);
		virtual void            assembleForce(const blue::ivector &f);
		virtual	void            assembleForce(const double d){}
		virtual void            assembleForceDof(const int dof, const double f);
		virtual blue::ivector         coordinates(const dofset::evaluation_time when=dofset::tna);
		virtual const blue::ivector	coordinates(const dofset::evaluation_time when=dofset::tna) const;
		blue::ivector&                displacement(const dofset::evaluation_time when=dofset::tna) { return U.displacement(when);}

		virtual	void            incrementSolution(const integrator& integ);

		virtual void            computeExplicitAcceleration();
		virtual void            computeExplicitVelocity();
        virtual void            computeExplicitEulerianAcceleration();
		
		void                    constrainDof(const int local);
		virtual blue::ivector&        force()						{return U.force();}
		virtual	void            forceReset()				{U.force().setZero();}
		virtual int             getID (const int local) const {return U.getID(local);}
		int                     getNDofs() const			{return U.ndof();}
		vectorDofset&           getUDS()					{return U;}
		const vectorDofset&     getUDS() const				{return U;}
		virtual int             getNConstrainedDofs() const	{return U.getNConstrainedDofs();}
		inline	nodetypeT       getNodetype() const			{return _Unode;}
		virtual bool            isDofConstrained(const int local) const {return U.isConstrained(local);}
		virtual bool            isOrphan()const				{return U.isUnused();}
		virtual void            localizeSolutionIncrement(const longvector& delta, partitionDofs pt=PARTITION_ALL);
		virtual void            printID(ostream &of=cout) const		{ of << U.getIDs();}
		virtual void            printDOFs(ostream &of=cout) const;
		virtual void            printForce(ostream &of=cout);
		virtual void            recoverFromBackupDofs();
		virtual void            setDeltaDofsToZero();
		virtual void            setDof(const int local, const double v);
		virtual void            setDofIncrement(const int local, const double v);
		virtual void            setInitialDof(const int local, const double v);
		virtual void            setReferenceCoordinates(const blue::ivector& c){refCoords = c; _currCoordinates = c+U.displacement();}
		blue::ivector&                velocity(const dofset::evaluation_time when=dofset::tna) { return U.velocity(when);}


		
		friend std::ostream&  operator<<(std::ostream &os, const Unode &nd);
		
	protected:
		vectorDofset	U;
		blue::ivector			_currCoordinates;
	};






class UPnode : public Unode
	{
		
	public:
		UPnode(const int label_, const blue::ivector& coor, const feliks::topology::cell *the0Cell=0);
		
		
		virtual	void			advanceIntegration(const integrator& integ, const double dt);
		virtual	void			incrementSolution(const integrator& integ);

		
		virtual void			assembleForce(const double f);
		virtual void			assembleForce(const blue::ivector& f){ Unode::assembleForce(f);}
		virtual void			assembleForceDof(const int dof, const double f);
		virtual void			constrainDof(const int local);

		
		virtual	void			forceReset()					{U.force().setZero(); P.force()=0.0;}
		virtual int				getID (const int local) const	{return local < 3 ? U.getID(local) : P.getID(local-3);}
		virtual int				getNConstrainedDofs() const		{return U.getNConstrainedDofs() + P.getNConstrainedDofs();}
		int						getNDofs() const				{return U.ndof() + P.ndof();}
		inline	nodetypeT		getNodetype() const				{return _UPnode;}
		scalarDofset&			getPDS()						{return P;}
		const scalarDofset&		getPDS()	const				{return P;}
		virtual bool			isDofConstrained(const int local) const {return local<3 ? U.isConstrained(local) : P.isConstrained(local);}
		virtual bool			isOrphan()const					{return U.isUnused() && P.isUnused();}
		virtual void			localizeSolutionIncrement(const longvector& delta,  partitionDofs pt=PARTITION_ALL);
		virtual void			printID(ostream &of=cout) const		{ of << U.getIDs() << P.getIDs();}
		virtual void			printDOFs(ostream &of=cout) const;
		virtual void			printForce(ostream &of=cout);
		void					recoverFromBackupDofs();
		virtual void			setDeltaDofsToZero();
		virtual void			setDof(const int local, const double v);
		virtual void			setDofIncrement(const int local, const double v);
		virtual void			setInitialDof(const int local, const double v);



		friend std::ostream&  operator<<(std::ostream &os, const UPnode &nd);

	protected:
		scalarDofset	P;
	};





class UDnode : public Unode
	{
	public:
		UDnode(const int label_, const blue::ivector& coor, const feliks::topology::cell *the0Cell=0);

		virtual	void			advanceIntegration(const integrator& integ, const double dt);
		virtual	void			incrementSolution(const integrator& integ);
		
		
		virtual void			assembleForce(const double f){}
		virtual void			assembleForce(const blue::ivector& f){ Unode::assembleForce(f);}
		virtual void			assembleForceDof(const int dof, const double f);
		virtual void			constrainDof(const int local);
		virtual	void			forceReset()					{U.force().setZero(); D.force()[0] = D.force()[1] = 0.0;}
		directorDofset&			getDDS()						{return D;}
		const directorDofset&	getDDS() const					{return D;}
		
		virtual int				getID (const int local) const	{return local < 3 ? U.getID(local) : D.getID(local-3);}
		virtual int				getNConstrainedDofs() const		{return U.getNConstrainedDofs() + D.getNConstrainedDofs();}
		int						getNDofs() const				{ return U.ndof() + D.ndof();}
		nodetypeT				getNodetype() const				{ return _UDnode;}
		virtual void			initializeRotationalDofs(const blue::ivector& director);
		virtual bool			isDofConstrained(const int local) const {return local<3 ? U.isConstrained(local) : D.isConstrained(local);}
		virtual bool			isOrphan()const					{return U.isUnused() && D.isUnused();}
		virtual void			localizeSolutionIncrement(const longvector& delta,  partitionDofs pt=PARTITION_ALL);
		virtual void			printID(ostream &of=cout) const		{ of << U.getIDs() << D.getIDs();}
		virtual void			printDOFs(ostream &of=cout) const;
		virtual void			printForce(ostream &of=cout);	
		void					recoverFromBackupDofs();
		virtual void			setDeltaDofsToZero();
		virtual void			setDofIncrement(const int local, const double v);
		virtual void			setInitialDof(const int local, const double v);



		friend std::ostream&  operator<<(std::ostream &os, const UDnode &nd);

	private:
		directorDofset	D;
		
	};



class Pnode : public node
	{
	public:
		Pnode(const int label_, const blue::ivector& coor, const feliks::topology::cell *the0Cell=0);

		virtual	void			advanceIntegration(const integrator& integ, const double dt);
		virtual	void			incrementSolution(const integrator& integ);

		virtual void			computeExplicitAcceleration();
		virtual void			computeExplicitVelocity();

        virtual void            computeExplicitEulerianAcceleration(){}; 

		virtual void			assembleForce(const double f);
		virtual	void			assembleForce(const blue::ivector& f){}
		virtual void			assembleForceDof(const int dof, const double f);
		virtual void			constrainDof(const int local);
		virtual	void			forceReset()					{P.force() = 0.0;}
		virtual int				getID (const int local) const	{return  P.getID(local);}
		scalarDofset&			getPDS()						{return P;}
		const scalarDofset&		getPDS() const					{return P;}
		virtual int				getNConstrainedDofs() const		{return P.getNConstrainedDofs();}
		int						getNDofs() const				{return P.ndof();}
		inline	nodetypeT		getNodetype() const				{return _Pnode;}
		virtual bool			isDofConstrained(const int local) const {return  P.isConstrained(local);}
		virtual bool			isOrphan()const					{return P.isUnused();}
		virtual void			localizeSolutionIncrement(const longvector& delta,  partitionDofs pt=PARTITION_ALL);
		virtual void			printDOFs(ostream &of=cout) const;
		virtual void			printID(ostream &of=cout) const		{ of << P.getIDs();}
		virtual void			printForce(ostream &of=cout);
		void					recoverFromBackupDofs();
		virtual void			setDeltaDofsToZero();
		virtual void			setDof(const double v);
		virtual void			setDofIncrement(const int local, const double v);
		virtual void			setInitialDof(const int local, const double v);


		friend std::ostream&  operator<<(std::ostream &os, const Pnode &nd); 

	
	protected:
		scalarDofset P;
	};



// a node holding a vector
template <int N>
class Vnode : public node
{
public:
	Vnode(const int label_, const blue::ivector& coor, const feliks::topology::cell *the0Cell=0);
	
	virtual	void			advanceIntegration(const integrator& integ, const double dt);
	virtual	void			incrementSolution(const integrator& integ);
	
	virtual void			computeExplicitAcceleration();
	virtual void			computeExplicitVelocity();
    virtual void            computeExplicitEulerianAcceleration(){}; 

	
	virtual void			assembleForce(const double* f);
	virtual void			assembleForceDof(const int dof, const double f);
	virtual void			constrainDof(const int local);
	virtual	void			forceReset()					{ds.resetForce();}
	virtual int				getID (const int local) const	{return ds .getID(local);}
	scalarDofset&			getVDS()						{return ds;}
	const scalarDofset&		getVDS() const					{return ds;}
	virtual int				getNConstrainedDofs() const		{return ds.getNConstrainedDofs();}
	int						getNDofs() const				{return ds.ndof();}
	inline	nodetypeT		getNodetype() const				{return _Vnode;}
	virtual bool			isDofConstrained(const int local) const {return  ds.isConstrained(local);}
	virtual bool			isOrphan()const					{return ds.isUnused();}
	virtual void			localizeSolutionIncrement(const longvector& delta,  partitionDofs pt=PARTITION_ALL);
	virtual void			printDOFs(ostream &of=cout) const;
	virtual void			printID(ostream &of=cout) const		{ of << ds.getIDs();}
	virtual void			printForce(ostream &of=cout);
	void					recoverFromBackupDofs();
	virtual void			setDeltaDofsToZero();
	virtual void			setDofIncrement(const int local, const double v);
	virtual void			setInitialDof(const int local, const double v);
	
	
	friend std::ostream&  operator<<(std::ostream &os, const Pnode &nd); 
	
	
protected:
	nvectorDofset<N> ds;
};




#endif
