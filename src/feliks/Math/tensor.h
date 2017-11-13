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
 * tensor.h
 *
 * Classes for tensors, vectors, and rotations in three dimensional Euclidean space.
 * operations on 3x3 matrices and 3x1 vectors
 *
 * ignacio romero, october 2002
 * revised november 2006, september 2008
 *
 *
 * The ivector class models vectors in R3. It supports the obvious access by the []
 * operator or by the () function call as well as the usual operators of sum,
 * difference, scalar multiplication and all the typical useful operations: norm, normalization,
 * triple product, angle with another vector,
 *
 */

#ifndef _blue_tensor_h
#define _blue_tensor_h

#include <cmath>
#include <vector>
#include <iostream>
#include <cassert>

namespace blue
{
    
    class iquaternion;
    class irotation;
    class itensor;
    
    
    // the ivector class has nothing to do with the stl vector. In fact, it is just a class
    // for vectors in R^3
    // the implementation of the functions is pretty straightforward, with no fancy memory
    // operations nor metaprogramming
    class ivector
    {
        
	public:
        ivector();
        ivector(const double a, const double b, const double c);
        ivector(const double *a);
        ivector(const ivector &v2);
        ivector(const std::vector<double>& v, const size_t start);
        
        
        ivector					operator-();
        inline ivector&			operator=(const  ivector &v);
        ivector&				operator+=(const ivector &v);
        ivector&				operator-=(const ivector &v);
        ivector&				operator*=(const double  a);
        inline double&			operator[](const int i);
        inline const double&	operator[](const int i) const;
        inline double&			operator()(const int i);
        inline const double&	operator()(const int i) const;
        
        friend std::ostream&	operator<<(std::ostream &os, const ivector& v);
        friend const ivector	operator+(const ivector &left, const ivector& right);
        friend const ivector	operator-(const ivector &left, const ivector& right);
        friend const ivector	operator*(const double  a,  const ivector& v);
        friend const ivector	operator*(const ivector& v, const double a);
        
        
        double					angleWith(const ivector& v2) const;
        double*					components();
        void                    changeSign();
        ivector					cross(const ivector &v) const;
        double					dot(const ivector &v) const;
        void					extractFrom(const iquaternion& q);
        void					extractFrom(const irotation& r);
        double					norm() const;
        double					squaredNorm() const;
        void					normalize();
        ivector					normalized();
        void					print(std::ostream &os=std::cout);
        void					scale(const double alpha);
        void					setRandom();
        void					setZero();
        static double			tripleProduct(const ivector &v1, const ivector &v2, const ivector &v3);
        
    private:
        double                  x[3];
        
    };
    
    
    // implementation of inline functions
    inline		double			ivector :: norm() const						{return sqrt( this->squaredNorm() );}
    inline		double			ivector :: squaredNorm() const				{return this->dot(*this);}
    inline		double*			ivector :: components()						{return x;}
    inline		void			ivector :: print(std::ostream &of)			{of << *this;}
    inline      double&			ivector :: operator[](const int i)			{assert(i<3); return x[i];}
    inline const double&		ivector :: operator[](const int i) const	{assert(i<3); return x[i];}
    inline		double&			ivector :: operator()(const int i)			{assert(i<3); return x[i];}
    inline const double&		ivector :: operator()(const int i) const	{assert(i<3); return x[i];}
    inline                      ivector :: ivector(const ivector& w)        {x[0]=w(0);	x[1]=w(1);x[2]=w(2);}
    inline ivector&             ivector :: operator=(const ivector &w)      {x[0]=w(0); x[1]=w(1);x[2]=w(2); return *this;}

    inline const ivector  operator+(const ivector &left, const ivector &right)
            { return ivector(left[0] + right[0] , left[1] + right[1] , left[2] + right[2] ); }
    
    inline const ivector  operator-(const ivector &left, const ivector &right)
            { return ivector(left[0] - right[0] , left[1] - right[1] , left[2] - right[2] ); }
        
    inline 	ivector&  ivector ::  operator+=(const ivector &v)
            { x[0] += v.x[0]; x[1] += v.x[1]; x[2] += v.x[2]; return *this;}
    
    
    inline const ivector operator*(const double a,  const ivector &w)
            { return ivector(a*w[0], a*w[1], a*w[2]); }
    
    inline const ivector operator*(const blue::ivector &w, const double a)
            { return ivector(a*w[0], a*w[1], a*w[2]); }
    
    

    // vectors in R2
    class vector2
    {
        
    public:
        vector2(){}
        vector2(const double x0, const double x1){ x[0]=x0, x[1] = x1;}
        vector2(const vector2& w);
        
        inline vector2&			operator=(const  vector2 &v);
        vector2&				operator+=(const vector2 &v);
        vector2&				operator-=(const vector2 &v);
        vector2&				operator*=(const double  a);
        inline double&			operator[](const int i)				{return x[i];}
        inline const double&	operator[](const int i) const		{return x[i];}
        inline double&			operator()(const int i)				{return x[i];}
        inline const double&	operator()(const int i) const		{return x[i];}
        
        friend std::ostream&	operator<<(std::ostream &os, const vector2 &v);
        friend const vector2	operator+(const vector2 &left, const vector2 &right);
        friend const vector2	operator-(const vector2 &left, const vector2 &right);
        friend const vector2	operator*(const double  a,  const vector2 &v);
        friend const vector2	operator*(const vector2 &v, const double a);
        
        inline void				print(std::ostream &of=std::cout);
        inline void				setZero(){x[0] = x[1] = 0.0;}
        
    private:
        double x[2];
    };
    
    
    inline vector2 :: vector2(const vector2& w)				{x[0] = w(0); x[1] = w(1);}
    inline void vector2 :: print(std::ostream &of) { of << x[0] << x[1];}
    
    
    inline vector2& vector2 ::operator=(const vector2 &w)	{x[0] = w(0); x[1] = w(1); return *this;}
    
    inline const vector2  operator+(const vector2 &left, const vector2 &right)
	{return vector2(left[0] + right[0] , left[1] + right[1]); }
    
    inline const vector2  operator-(const vector2 &left, const vector2 &right)
	{return vector2(left[0] - right[0] , left[1] - right[1]); }
    
    inline const vector2  operator*(const vector2 &v, const double a)
	{return vector2(v[0]*a , v[1]*a); }
    
    inline const vector2  operator*(const double a, const vector2 &v)
	{return vector2(v[0]*a , v[1]*a); }
    
    inline 	vector2&  vector2 ::  operator+=(const vector2 &v)
	{ x[0] += v.x[0]; x[1] += v.x[1]; return *this;}
    
    
    
    
    
    
    // General rank-2 tensor for 2 or three dimensional problems. Tensors on reals
    //
    class itensor
    {
        
    public:
        itensor();
        itensor(const itensor &t);
        itensor(const ivector& col1, const ivector& col2, const ivector& col3);
        itensor(const double a00, const double a01, const double a02,
                const double a10, const double a11, const double a12,
                const double a20, const double a21, const double a22);
        
        
        virtual ~itensor(){};
        
        static const itensor  dyadic(const ivector&a, const ivector& b);
        static const itensor  identity();
        
        
        itensor&        operator=( const itensor &t);
        itensor 		operator-();
        itensor&        operator+=(const itensor &t);
        itensor&        operator-=(const itensor &t);
        itensor&        operator*=(const double  alpha);
        double&         operator()(const int i, const int j)	 {return a[i][j];}
        const double&   operator()(const int i, const int j)const{return a[i][j];}
        
        friend itensor  operator+(const itensor &t1, const itensor &t2);
        friend itensor  operator-(const itensor &t1, const itensor &t2);
        friend itensor  operator*(const itensor &t1, const itensor &t2);
        friend itensor  operator*(const itensor &t , const double   a);
        friend itensor  operator*(const double  a  , const itensor &t);
        friend ivector  operator*(const itensor &t , const ivector &v);
        
        friend std::ostream& operator<<(std::ostream &os, const itensor &t);
        
        void            addDyadic(const ivector &a, const ivector &b);
        virtual void    addSymmetrizedDyadic(const ivector &a, const ivector &b);
        ivector			axialVector() const;
        void			beDexp(const ivector& a);
        void			beDexpinv(const ivector &a);
        void            beSkew(const ivector &v);
        void            beSkew();
        virtual ivector col(const int n) const;
        virtual double  contract(const itensor &t) const;
        virtual double  determinant() const;
        virtual double  dot(const itensor &t) const;
        virtual double  invariant1() const;
        virtual double  invariant2() const;
        virtual double  invariant3() const;
        itensor			inverse() const;
        virtual double  invert();
        virtual double  norm() const;
        virtual ivector row(const int n) const;
        virtual void	setRandom();
        virtual void	setZero();
        static  itensor skew(const ivector& v);
        double			squaredNorm() const;
        virtual double  trace() const;
        virtual itensor transpose() const;
        
        virtual	void jacobi_rot(const double s, const double tau,  const int i, const  int j, const  int k, const int l);
        
        
    protected:
        
        double a[3][3];
        
    };
    
    
    // implementation of inline functions
    inline  double	itensor :: squaredNorm() const			{return contract(*this);}
    
    
    
    class istensor : public itensor
    {
        
        
    public:
        istensor();
        istensor(const istensor &T);
        
        static const istensor  identity();
        static const istensor  deviatoricPart(const istensor& t);
        static const istensor  random();
        static const istensor  symmetricPartOf(const itensor& t);
        static const istensor  tensorTimesTensorTransposed(const itensor& t);
        static const istensor  tensorTransposedTimesTensor(const itensor& t);
        static const istensor  FtCF(const itensor& F, const istensor& C);
        static const istensor  FSFt(const itensor& F, const istensor& S);
        
		
        istensor&        operator=(const  istensor &t);
        istensor 		 operator-();
        istensor&        operator+=(const istensor &t);
        istensor&        operator-=(const istensor &t);
        istensor&        operator*=(const double alpha);
        
        friend istensor  operator+(const istensor &t1, const istensor &t2);
        friend itensor   operator+(const istensor &t1, const itensor &t2);
        friend itensor   operator+(const itensor  &t1, const istensor &t2);
        friend istensor  operator-(const istensor &t1, const istensor &t2);
        friend istensor  operator*(const istensor &t1, const istensor &t2);
        friend itensor   operator*(const itensor  &t1, const istensor &t2);
        friend itensor   operator*(const istensor &t1, const itensor  &t2);
        friend istensor  operator*(const istensor &t , const double   a);
        friend istensor  operator*(const double   a  , const istensor &t);
        friend ivector   operator*(const istensor &t , const ivector &v);
        
        void            addScaledVdyadicV(const double alpha, const ivector& V);
        void            bePushForwardContraContra(const istensor &S, const itensor &F);
        void            beSymmetricPartOf(const itensor &t);
        const ivector   eigenvalues() const;
        const istensor	inverse() const;
        double          maxEigenvalue() const;
        bool            maxPrincipalDirection(ivector& v) const;
        void            setRandom();
        void            spectralDecomposition(ivector ivectors[3], ivector &evalues) const;
        
    private:
        istensor&        operator=(const  itensor  &t);
        istensor(const itensor  &T);
    };
    
    
    
    class skewtensor : public itensor
    {
		
    private:
        skewtensor&        operator=(const  itensor  &t);
        skewtensor(const itensor  &T);
    };
    
    
    
    
    
    class iquaternion
    {
    public:
        
        iquaternion();
        iquaternion(const double q0, const double q1, const double q2, const double q3);
        iquaternion(const double* q);
        iquaternion(const iquaternion& q_);
        iquaternion(const irotation &m);
        iquaternion(const ivector &v);
        iquaternion& operator=(const iquaternion &rhs);
        iquaternion& operator=(const irotation &rhs);
        
        
        double*             components(){return q;}
        iquaternion         conjugate();
        void                extractFromRotation(const irotation &m);
        double              norm() const;
        void                normalize();
        void                toIdentity();
        
        friend ivector		operator*(const iquaternion &quat , const ivector &v);
        friend iquaternion	operator*(const iquaternion &q1, const iquaternion &q2);
        friend std::ostream& operator<<(std::ostream &os, const iquaternion &t);
        
        double&             x();
        double&             y();
        double&             z();
        double&             w();
        const double&       x() const;
        const double&       y() const;
        const double&       z() const;
        const double&       w() const;
        
        
    private:
        double q[4];
        friend class irotation;
    };
    
    
    
    
    class irotation : public itensor
    {
        
    public:
        irotation();
        irotation(const ivector &theta);
        irotation(const iquaternion &q);
        irotation(const double a, const double b, const double c);
        irotation(const itensor& t);
        
        ivector		rotationVector() const;
        void		beRotationWithoutDrill(const ivector &v);
        void		beRotationWithoutDrill(const ivector &from, const ivector& to);
        void		computeFrom(const iquaternion &q);
        
    
    private:
        
        //prevent some functions by declaring them private and not defining them
        friend irotation   operator+(const irotation &t1, const irotation &t2);
        friend irotation   operator-(const irotation &t1, const irotation &t2);
        friend irotation   operator*(const irotation &t1, const double   a);
        friend irotation   operator*(const double  a    , const irotation &t);
        
        irotation&  operator*=(const double alpha);
    
    };
    
}

#endif
