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
/* tensor.cpp
 * i. romero
 *
 * provides implementation for the most typical operations involving vectors, tensors, symmetric tensors,
 * rotation tensors, and iquaternions.
 *
 * The resulting objects can be manipulated with ease and written in a very compact form, almost as
 * one would write it in math. For certain operations (e.g., A = B*C*D, all being tensors) there
 * will be a significant overhead when compared with all loop implementations due to temporaries.
 * There are ways to eliminate this (using Expression Templates as in uBLAS, Boost, Blitz, ...) but always at the expense
 * of some syntax modification. In the future we should look into this more carefully, however,
 * for the moment, we will keep on working with these functions.
 *
 * Instead of using loops, often we have chosen to write operations explicitly. This is because for very short
 * vectors and tensors, a great deal of time is wasted in the comparison, incrementation, ..., relatively
 * to the time spent in the math operations. Moreover, writing it out explicitly helps some processors to
 * vectorize the operations.
 *
 * We could have used a templated definition for vectors of dimension 1, 2, or 3. Instead, we have simply,
 * for the moment, implemented vectors of fixed size 3, making sure the third coordinate is 0 always, in
 * case ndm = 2
 */

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <utility>

#include <cmath>
#include <cstdlib>
#include <vector>
#include "Math/tensor.h"
#include "Math/statistics.h"

using namespace blue;



ivector :: ivector()
{
    x[0] = x[1] = x[2] = 0.0;
}




ivector :: ivector(const double* a)
{
    x[0] = a[0];
    x[1] = a[1];
    x[2] = a[2];
}




ivector :: ivector(const double alpha, const double beta, const double gamma)
{
	x[0] = alpha;
	x[1] = beta;
	x[2] = gamma;
}




ivector :: ivector(const std::vector<double>& v, const size_t start)
{
	x[0] = v[start];
	x[1] = v[start+1];
	x[2] = v[start+2];
}


ivector ivector ::	operator-()
{
	return ivector(-x[0], -x[1], -x[2]);
}



ivector& ivector ::  operator-=(const ivector &w)
{
	x[0] -= w(0);
	x[1] -= w(1);
	x[2] -= w(2);
	
	return *this;
}



ivector&  ivector :: operator*=(const double a)
{
	x[0] *= a;
	x[1] *= a;
	x[2] *= a;
	return *this;
}



namespace blue
{
    
    std::ostream& operator<<(std::ostream &os, const ivector &v)
    {
        os  << std::setw(12) << std::setprecision(8)
        << "[ " << v.x[0] << " , " << v.x[1] << " , " << v.x[2]
        << " ]" << std::flush;
        
        return os;
    }
    
}




double ivector :: angleWith(const ivector& v2) const
{
	return ( dot(v2) / ( norm() * v2.norm() ) );
}




void ivector :: changeSign()
{
    x[0] = -x[0];
    x[1] = -x[1];
    x[2] = -x[2];
}




ivector ivector :: cross(const ivector &w) const
{
	return ivector(x[1]*w[2] - x[2]*w[1],
				   x[2]*w[0] - x[0]*w[2],
				   x[0]*w[1] - x[1]*w[0]);
}


double ivector :: dot(const ivector &w) const
{
	return x[0]*w.x[0] + x[1]*w.x[1] + x[2]*w.x[2];
}



void ivector :: extractFrom(const iquaternion& quat)
{
	double qnorm  = sqrt(quat.x()*quat.x() + quat.y()*quat.y() + quat.z()*quat.z());
    double mm     = (1.0 < qnorm) ? 1.0 : qnorm;
	double rotnr2 = asin( mm );
	double rotnrm = rotnr2 * 2.0;
	
	double rotfac = (qnorm > 1.0e-10) ? rotnrm / qnorm : 2.0;
	
    x[0] = rotfac * quat.x();
    x[1] = rotfac * quat.y();
    x[2] = rotfac * quat.z();
}



void ivector ::	extractFrom(const irotation& r)
{
	extractFrom( iquaternion(r) );
}



void ivector ::  normalize()
{
	(*this) *= (1.0/norm());
}




ivector ivector ::  normalized()
{
	double in(1.0/norm());
	return ivector( in * (*this) );
}




void ivector :: scale(const double alpha)
{
	x[0] *= alpha;
	x[1] *= alpha;
	x[2] *= alpha;
}




void ivector :: setRandom()
{
	x[0] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
	x[1] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
	x[2] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
}



void ivector :: setZero()
{
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
}



double ivector :: tripleProduct(const ivector& v1, const ivector &v2, const ivector &v3)
{
	const ivector tmp(v1.cross(v2));
	return v3.dot(tmp);
}




itensor :: itensor()
{
	a[0][0] = a[0][1] = a[0][2] = 0.0;
	a[1][0] = a[1][1] = a[1][2] = 0.0;
	a[2][0] = a[2][1] = a[2][2] = 0.0;
}



itensor :: itensor(const itensor &t)
{
	a[0][0] = t.a[0][0]; a[0][1] = t.a[0][1]; a[0][2] = t.a[0][2];
	a[1][0] = t.a[1][0]; a[1][1] = t.a[1][1]; a[1][2] = t.a[1][2];
	a[2][0] = t.a[2][0]; a[2][1] = t.a[2][1]; a[2][2] = t.a[2][2];
}


itensor :: itensor(const ivector& col1, const ivector& col2, const ivector& col3)
{
	const int ndm(3);
	for (int i=0; i<ndm; i++)
	{
		a[i][0] = col1[i];
		a[i][1] = col2[i];
		a[i][2] = col3[i];
	}
}




itensor :: itensor(const double a00, const double a01, const double a02,
				   const double a10, const double a11, const double a12,
				   const double a20, const double a21, const double a22)
{
	a[0][0] = a00; 	a[0][1] = a01; 	a[0][2] = a02;
	a[1][0] = a10; 	a[1][1] = a11; 	a[1][2] = a12;
	a[2][0] = a20; 	a[2][1] = a21; 	a[2][2] = a22;
}




void  itensor :: addDyadic(const ivector &m, const ivector &n)
{
    a[0][0] += m[0]*n[0];
    a[0][1] += m[0]*n[1];
    a[0][2] += m[0]*n[2];
    a[1][0] += m[1]*n[0];
    a[1][1] += m[1]*n[1];
    a[1][2] += m[1]*n[2];
    a[2][0] += m[2]*n[0];
    a[2][1] += m[2]*n[1];
    a[2][2] += m[2]*n[2];
    /*
     const int ndm(3);
     for (int i=0; i<ndm; i++)
     for (int j=0; j<ndm; j++)
     a[i][j] += m[i]*n[j];
     */
}





void itensor :: addSymmetrizedDyadic(const ivector &m, const ivector &n)
{
	double tmp;
	
	a[0][0] +=  m[0]*n[0];
	a[1][1] +=  m[1]*n[1];
	
	tmp      =  0.5*(m[0]*n[1]+m[1]*n[0]);
	a[0][1] +=  tmp;
	a[1][0] +=  tmp;
	
	tmp      = 0.5*(m[0]*n[2]+m[2]*n[0]);
	a[0][2] += tmp;
	a[2][0] += tmp;
    
	tmp      = 0.5*(m[1]*n[2]+m[2]*n[1]);
	a[1][2] += tmp;
	a[2][1] += tmp;
    
	a[2][2] += m[2]*n[2];
}



ivector itensor ::	axialVector() const
{
	itensor s = 0.5*(*this - this->transpose());
	return ivector( -s(1,2), -s(2,0), -s(0,1) );
}


/*
 void itensor ::	beAdjoint()
 {
 itensor A( *this );
 
 a[0][0] =  A(1,1)*A(2,2)-A(1,2)*A(2,1);
 a[1][0] = -A(2,1)*A(0,2)+A(2,2)*A(0,1);
 a[2][0] =  A(0,1)*A(1,2)-A(0,2)*A(1,1);
 
 a[0][1] = -A(1,2)*A(2,0)+A(1,0)*A(2,2);
 a[1][1] =  A(2,2)*A(0,0)-A(2,0)*A(0,2);
 a[2][1] = -A(0,2)*A(1,0)+A(0,0)*A(1,2);
 
 a[0][2] =  A(1,0)*A(2,1)-A(1,1)*A(2,0);
 a[1][2] = -A(2,0)*A(0,1)+A(2,1)*A(0,0);
 a[2][2] =  A(0,0)*A(1,1)-A(0,1)*A(1,0);
 }
 */




// linearized exponential map dexp(psi)
/*
 *  Notes:      d = norm(psi)
 *              n = psi/d
 *              nn= n otimes n
 *
 *              h = nn + sin(d)/d (eye-nn) + (1-cos(d))/d hatn,
 */
void itensor :: beDexp(const ivector &psi)
{
    double  s, c, d, z;
	
    double b = psi.dot(psi);
	double q;
    
    // limit case
    if (b < 1.0e-12)
    {
        s = (1.0 - 0.05*b)/6.0;
        c =  0.5 - b*(1.0 - b/30.0)/24.0;
        q =  1.0 - b*s;
    }
    else
    {
        d = sqrt(b);
        z = sin(d)/d;
        c = (1.0 - cos(d))/b;
        s = (1.0-z)/b;
        q = z;
    }
	
	itensor hn;
    hn.beSkew(psi);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
            a[i][j] = s*psi[i]*psi[j] + c*hn(i,j);
        a[i][i] +=  q;
    }
}






// inverse of the linearized exp operator
void itensor :: beDexpinv(const ivector &theta)
{
	int i   , j;
	double nn, norm;
	double f, c;
	itensor htheta;
	
	nn = theta.dot(theta);
	
	if (nn < 1e-12)
	{
		// new expansion
		f = 1.0 - nn/12.0;
		c = 1.0/12.0 + nn/720.0;
	}
	else
	{
		norm = sqrt(nn);
		f    = 0.5*norm/(tan(0.5*norm));
		c    = (1.0 - f)/nn;
	}
	
	htheta.beSkew(theta);
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
			a[i][j] = c*theta[i]*theta[j] - 0.5*htheta(i,j);
		a[i][i] += f;
	}
}


void itensor :: beSkew(const ivector &v)
{
	a[0][0] =   0.0;
	a[1][1] =   0.0;
	a[2][2] =   0.0;
	
    a[0][1] = -v[2];
    a[1][0] =  v[2];
	
    a[0][2] =  v[1];
	a[2][0] = -v[1];
    
    a[1][2] = -v[0];
    a[2][1] =  v[0];
}




/* test code to convert a second order generic tensor to
 * a antisymmetric one */
void itensor :: beSkew()
{
	a[0][0] =   0.0;
	a[1][1] =   0.0;
	a[2][2] =   0.0;
	
    a[0][1] = 0.5*(a[0][1]-a[1][0]);
    a[0][2] = 0.5*(a[0][2]-a[2][0]);
    a[1][2] = 0.5*(a[1][2]-a[2][1]);
	
    a[1][0] = -a[0][1];
    a[2][0] = -a[0][2];
    a[2][1] = -a[1][2];
}



ivector itensor :: col(const int n) const
{
	return ivector( a[0][n], a[1][n], a[2][n] );
}


double itensor :: contract(const itensor &U) const
{
	double d=0.0;
	const int ndm(3);
    
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			d += a[i][j] * U.a[i][j];
	return d;
}




double itensor :: determinant() const
{
	double d;
	
	d = +   a[0][0]*a[1][1]*a[2][2]
    +   a[1][0]*a[2][1]*a[0][2]
    +   a[2][0]*a[0][1]*a[1][2]
    -   a[0][2]*a[1][1]*a[2][0]
    -   a[0][1]*a[1][0]*a[2][2]
    -   a[0][0]*a[2][1]*a[1][2];
    
	return d;
}



double itensor :: dot(const itensor &U) const
{
	return contract(U);
}




const itensor  itensor :: dyadic(const ivector&a, const ivector& b)
{
	return itensor (a(0)*b(0), a(0)*b(1), a(0)*b(2),
					a(1)*b(0), a(1)*b(1), a(1)*b(2),
					a(2)*b(0), a(2)*b(1), a(2)*b(2));
}



const itensor  itensor :: identity()
{
	return itensor( 1.0, 0.0, 0.0,
                   0.0, 1.0, 0.0,
                   0.0, 0.0, 1.0);
}


double itensor :: invariant1() const
{
	return trace();
}


double itensor :: invariant2() const
{
	itensor tt( (*this)*(*this));
	double  z( trace() );
	return 0.5*( z*z - tt.trace());
}



double itensor :: invariant3() const
{
	return determinant();
}


itensor itensor :: inverse() const
{
	double det = determinant();
	itensor m;
	
	if (det == 0.0)
		std::cout << "\n ERROR in tensor inverse. Tensor is singular" << std::endl;
	else
    {
		double idet = 1.0/det;
        
        m.a[0][0] = idet*(a[1][1]*a[2][2] - a[2][1]*a[1][2]);
        m.a[0][1] = idet*(a[0][2]*a[2][1] - a[0][1]*a[2][2]);
        m.a[0][2] = idet*(a[0][1]*a[1][2] - a[0][2]*a[1][1]);
      	
        m.a[1][0] = idet*(a[1][2]*a[2][0] - a[1][0]*a[2][2]);
        m.a[1][1] = idet*(a[0][0]*a[2][2] - a[0][2]*a[2][0]);
        m.a[1][2] = idet*(a[1][0]*a[0][2] - a[0][0]*a[1][2]);
		
        m.a[2][0] = idet*(a[1][0]*a[2][1] - a[1][1]*a[2][0]);
        m.a[2][1] = idet*(a[0][1]*a[2][0] - a[0][0]*a[2][1]);
        m.a[2][2] = idet*(a[0][0]*a[1][1] - a[0][1]*a[1][0]);
	}
	return m;
}


double itensor :: invert()
{
	double det = determinant();
	itensor m;
	
	if (det == 0.0) std::cout << "\n ERROR in tensor inverse. Tensor is singular" << std::endl;
	else
    {
		double idet = 1.0/det;
		
        m.a[0][0] = idet*(a[1][1]*a[2][2] - a[2][1]*a[1][2]);
        m.a[0][1] = idet*(a[0][2]*a[2][1] - a[0][1]*a[2][2]);
        m.a[0][2] = idet*(a[0][1]*a[1][2] - a[0][2]*a[1][1]);
        
        m.a[1][0] = idet*(a[1][2]*a[2][0] - a[1][0]*a[2][2]);
        m.a[1][1] = idet*(a[0][0]*a[2][2] - a[0][2]*a[2][0]);
        m.a[1][2] = idet*(a[1][0]*a[0][2] - a[0][0]*a[1][2]);
        
        m.a[2][0] = idet*(a[1][0]*a[2][1] - a[1][1]*a[2][0]);
        m.a[2][1] = idet*(a[0][1]*a[2][0] - a[0][0]*a[2][1]);
        m.a[2][2] = idet*(a[0][0]*a[1][1] - a[0][1]*a[1][0]);
	}
	*this = m;
	
	return det;
}



double itensor :: norm() const
{
	return sqrt( contract(*this));
}



ivector itensor :: row(const int n) const
{
	return ivector( a[n][0], a[n][1], a[n][2] );
}




void itensor :: setZero()
{
	const int ndm(3);
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] = 0.0;
}



itensor& itensor :: operator=(const itensor &t)
{
	const int ndm(3);
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] = t.a[i][j];
    
	return *this;
}



itensor itensor :: operator-()
{
	const int ndm(3);
	itensor ret;
	
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			ret(i,j) = -a[i][j];
	return ret;
}



itensor&  itensor ::  operator+=(const itensor &t)
{
	const int ndm(3);
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] += t.a[i][j];
	return *this;
}


itensor&  itensor ::  operator-=(const itensor &t)
{
	const int ndm(3);
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] -= t.a[i][j];
	return *this;
}




itensor&  itensor ::  operator*=(const double alpha)
{
	const int ndm(3);
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] *= alpha;
	return *this;
}


namespace blue
{
    
    std::ostream& operator<<(std::ostream &os, const itensor &t)
    {
        os  << std::setw(12) << std::setprecision(8)
        << "[ " << t.a[0][0] << " , " << t.a[0][1]  << " , " << t.a[0][2] << " ]\n";
        
        os  << std::setw(12) << std::setprecision(8)
        << "[ " << t.a[1][0] << " , " << t.a[1][1]  << " , " << t.a[1][2] << " ]\n";
        
        os  << std::setw(12) << std::setprecision(8)
        << "[ " << t.a[2][0] << " , " << t.a[2][1] 	<< " , " << t.a[2][2] << " ]\n";
        
        return os;
    }
    
}




itensor itensor :: skew(const ivector& v)
{
	itensor vhat; vhat.setZero();
	vhat(0,1) = -v(2);
	vhat(0,2) =  v(1);
	vhat(1,0) =  v(2);
	vhat(1,2) = -v(0);
	vhat(2,0) = -v(1);
	vhat(2,1) =  v(0);
	return vhat;
}




void itensor :: setRandom()
{
	const int ndm(3);
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] = random::uniform(-1.0,1.0);
}




double itensor :: trace() const
{
	double t=0.0;
	const int ndm(3);
    
	for (int i=0; i<ndm; i++) t += a[i][i];
	return t;
}



itensor itensor :: transpose() const
{
	itensor m;
	const int ndm(3);
    
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			m.a[i][j] = a[j][i];
	
	return m;
}


namespace blue
{
    
    
    itensor operator+(const itensor &t1, const itensor &t2)
    {
        itensor ret(t1);
        ret += t2;
        return ret;
    }
    
    
    itensor operator-(const itensor &t1, const itensor &t2)
    {
        itensor ret(t1);
        ret -= t2;
        return ret;
    }
    
    
    
    itensor  operator*(const itensor &T, const itensor &U)
    {
        itensor ret;
        const int ndm(3);
        
        for (int i=0; i<ndm; i++)
            for (int j=0; j<ndm; j++)
                for (int k=0; k<ndm; k++)
                    ret.a[i][j] += T.a[i][k] * U.a[k][j];
        
        return ret;
    }
    
    
    ivector  operator*(const itensor &t , const ivector &v)
    {
        ivector w;
        const int ndm(3);
        
        for (int i=0; i<ndm; i++)
            for (int j=0; j<ndm; j++)
                w[i] += t.a[i][j] * v[j];
        return w;
    }
    
    
    itensor  operator*(const itensor &t , const double alpha)
    {
        itensor ret(t);
        ret *= alpha;
        return ret;
    }
    
    
    itensor  operator*(const double alpha, const itensor &t)
    {
        return t*alpha;
    }
    
}

/* implementation of symmetric, second order tensor. This implementation should be
 changed in the future, to operate only with the upper or lower part of the tensor. For the
 moment, I just copy the itensor implementation and make it work
 */

istensor :: istensor() : itensor()
{
}



// for this constructor, we assume that T is symmetric
istensor :: istensor(const istensor &T)
{
	const int ndm(3);
    
	for (int i=0; i<ndm; i++)
	{
		a[i][i] = T(i,i);
		for (int j=i+1; j<ndm; j++)
			a[i][j] = a[j][i] = T(i,j);
	}
}


// for this constructor, we assume that T is symmetric
istensor :: istensor(const itensor &T)
{
	*this = symmetricPartOf(T);
}



void istensor :: addScaledVdyadicV(const double alpha, const ivector& V)
{
    const int ndm(3);
	
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
            a[i][j] += alpha * V(i)*V(j);
}


void istensor :: bePushForwardContraContra(const istensor &S, const itensor &F)
{
	const int ndm(3);
	double R[3][3];
	
	// multiply F*S using symmetry of S
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
		{
			R[i][j] = 0.0;
			for (int k=0; k<ndm; k++)
				R[i][j] += F(i,k) * S(j,k);
		}
    
	// multiply R*F^T
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
        {
            a[i][j] = 0.0;
            for (int k=0; k<ndm; k++)
                a[i][j] += R[i][k] * F(j,k);
        }
}



void istensor :: beSymmetricPartOf(const itensor &t)
{
	double tmp;
	
	a[0][0]  =  t(0,0);
	a[1][1]  =  t(1,1);
	
	tmp      =  0.5*(t(0,1)+t(1,0));
	a[0][1]  =  tmp;
	a[1][0]  =  tmp;
	
    tmp      = 0.5*(t(0,2)+t(2,0));
    a[0][2]  = tmp;
    a[2][0]  = tmp;
	
    tmp      = 0.5*(t(1,2)+t(2,1));
    a[1][2]  = tmp;
    a[2][1]  = tmp;
	
    a[2][2]  = t(2,2);
}



const istensor istensor :: identity()
{
    istensor t;
    t.setZero();
    t(0,0) = t(1,1) = t(2,2) = 1.0;
    return t;
}



const istensor istensor :: inverse() const
{
    double det = determinant();
	istensor m;
	
	if (det == 0.0)
		std::cout << "\n ERROR in tensor inverse. Tensor is singular" << std::endl;
	else
    {
		double idet = 1.0/det;
        
        m.a[0][0] = idet*(a[1][1]*a[2][2] - a[2][1]*a[1][2]);
        m.a[0][1] = idet*(a[0][2]*a[2][1] - a[0][1]*a[2][2]);
        m.a[0][2] = idet*(a[0][1]*a[1][2] - a[0][2]*a[1][1]);
      	
        m.a[1][1] = idet*(a[0][0]*a[2][2] - a[0][2]*a[2][0]);
        m.a[1][2] = idet*(a[1][0]*a[0][2] - a[0][0]*a[1][2]);
		
        m.a[2][2] = idet*(a[0][0]*a[1][1] - a[0][1]*a[1][0]);
        
        m(1,0) = m(0,1);
        m(2,0) = m(0,2);
        m(2,1) = m(1,2);
	}
    
	return m;
}


// compute the three eigenvalues, and take the maximum.
// Since the matrix is symmetric, the three roots of the 3rd order
// polynomial are real
double istensor :: maxEigenvalue() const
{
    ivector eval = (*this).eigenvalues();
    return eval(2);
}


//eigenvalues are sorted from smallest to largest
const ivector istensor :: eigenvalues() const
{
    ivector eval;
    
    
	// polynomial is x^3 + a x^2  + b  x + c  = 0
	//              -l^3 + I1 x^2 - I2 x + I3 = 0
	double a, b, c;
	a = -invariant1();
	b =  invariant2();
	c = -invariant3();
	
	if ( std::abs(a) + std::abs(b) + std::abs(c) <= 1e-8)
		eval.setZero();
    
	else
	{
		double  pi  = 3.14159265358;
		double	oot = 1.0/3.0;
		double  opf = 1.5;
		
		double p = (3.0*b-a*a)*oot;
		double q = c + 2.0*a*a*a/27.0-a*b*oot;
		
		// three real roots
		double  ap     = std::abs(p);
		double	cosphi = -0.5*q/std::pow(ap*oot, opf);
		double  phi;
		if (cosphi >= 1.0)
			phi = 0.0;
		else if (cosphi <= -1.0)
			phi = pi;
		else
			phi = std::acos(cosphi);
		
		if ( std::isnan(cosphi) )
			phi = pi;
		
		double	cf = 2.0*std::sqrt(ap*oot);
        
        double l0, l1, l2;
		l0 = -a*oot + cf*std::cos(     phi*oot);
		l1 = -a*oot - cf*std::cos((phi-pi)*oot);
		l2 = -a*oot - cf*std::cos((phi+pi)*oot);
		
        if (l0>l1) std::swap(l0, l1);
        if (l1>l2) std::swap(l1, l2);
        if (l0>l1) std::swap(l0, l1);
        
        eval(0) = l0;
        eval(1) = l1;
        eval(2) = l2;
	}
	
	return eval;
}



// copied from Sancho and Planas, which I believed copied it from feap
bool istensor :: maxPrincipalDirection(ivector& v) const
{
	istensor copy(*this);
	itensor  evec;
	double   eval[3];
	int nrot;
	bool ret(true);
	
	int i,j,ip,iq;
	double tresh,theta,tau,t,sm,s,h,g,c;
	
	int n = 3;
	double b[3], z[3];
	
	evec = itensor::identity();
	
	for (ip=0; ip<n; ip++)
	{
		b[ip] = eval[ip] = copy.a[ip][ip];
		z[ip] = 0.0;
	}
	
	nrot=0;
	for (i=1;i<=50;i++)
	{
		sm=0.0;
		for (ip=0;ip<n-1;ip++)
		{
			for (iq=ip+1;iq<n;iq++)
				sm += fabs(copy.a[ip][iq]);
		}
        
		if (sm == 0.0)
			break;
		
		if (i < 4)
			tresh = 0.2*sm/(n*n);
		
		else
			tresh = 0.0;
		
		for (ip=0; ip<n-1; ip++)
		{
			for (iq=ip+1; iq<n; iq++)
			{
				g = 100.0*fabs(copy.a[ip][iq]);
				if (i > 4 && (fabs(eval[ip])+g) == fabs(eval[ip]) && (fabs(eval[iq])+g) == fabs(eval[iq]))
					copy.a[ip][iq] = 0.0;
				else if (fabs(copy.a[ip][iq]) > tresh)
				{
					h=eval[iq]-eval[ip];
					if ((fabs(h)+g) == fabs(h))
						t=(copy.a[ip][iq])/h;
					else
					{
						theta=0.5*h/(copy.a[ip][iq]);
						t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if (theta < 0.0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*copy.a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					eval[ip] -= h;
					eval[iq] += h;
					copy.a[ip][iq]=0.0;
					for (j=0; j<ip; j++) 			copy.jacobi_rot(s, tau, j, ip, j, iq);
					for (j=ip+1; j<iq; j++)			copy.jacobi_rot(s, tau, ip, j, j, iq);
					for (j=iq+1; j<n; j++) 			copy.jacobi_rot(s, tau, ip, j, iq, j);
					for (j=0; j<n; j++)				evec.jacobi_rot(s, tau, j, ip, j, iq);
					++nrot;
				}
			}
		}
		
		for (ip=0;ip<n;ip++)
		{
			b[ip]    += z[ip];
			eval[ip]  = b[ip];
			z[ip]     = 0.0;
		}
	}
	if (i == 50) ret = false;
	
	
	// extract the maximum direction by comparing the three eigenvalues
	double max(eval[0]);
	i = 0;
	
	if ( eval[1] > eval[0] ) {max = eval[1]; i=1;}
	if ( eval[2] > max     ) {i=2;}
	
	v(0) = evec(0, i);
	v(1) = evec(1, i);
	v(2) = evec(2, i);
	
	return ret;
}


void itensor :: jacobi_rot(const double s, const double tau,  const int i, const  int j, const  int k, const int l)
{
	double g,h;
	
	g = a[i][j];
	h = a[k][l];
	a[i][j] = g - s*(h+g*tau);
	a[k][l] = h + s*(g-h*tau);
}




/* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
 domain Java Matrix library JAMA. */
#include <math.h>

#ifdef MAX
#undef MAX
#endif

#define MAX(a, b) ((a)>(b)?(a):(b))


static double hypot2(double x, double y) {
	return sqrt(x*x+y*y);
}


// Symmetric Householder reduction to tridiagonal form.
static void tred2(double V[3][3], double d[3], double e[3])
{
	const int n=3;
	
	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.
	
	for (int j = 0; j < n; j++)
    {
		d[j] = V[n-1][j];
	}
	
	// Householder reduction to tridiagonal form.
	
	for (int i = n-1; i > 0; i--)
    {
		
		// Scale to avoid under/overflow.
		
		double scale = 0.0;
		double h = 0.0;
		for (int k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i-1];
			for (int j = 0; j < i; j++) {
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		} else {
			
			// Generate Householder vector.
			
			for (int k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i-1];
			double g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i-1] = f - g;
			for (int j = 0; j < i; j++) {
				e[j] = 0.0;
			}
			
			// Apply similarity transformation to remaining columns.
			
			for (int j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j+1; k <= i-1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for (int j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for (int j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for (int j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i-1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i-1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}
	
	// Accumulate transformations.
	
	for (int i = 0; i < n-1; i++) {
		V[n-1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i+1];
		if (h != 0.0) {
			for (int k = 0; k <= i; k++) {
				d[k] = V[k][i+1] / h;
			}
			for (int j = 0; j <= i; j++) {
				double g = 0.0;
				for (int k = 0; k <= i; k++) {
					g += V[k][i+1] * V[k][j];
				}
				for (int k = 0; k <= i; k++) {
					V[k][j] -= g * d[k];
				}
			}
		}
		for (int k = 0; k <= i; k++) {
			V[k][i+1] = 0.0;
		}
	}
	for (int j = 0; j < n; j++) {
		d[j] = V[n-1][j];
		V[n-1][j] = 0.0;
	}
	V[n-1][n-1] = 1.0;
	e[0] = 0.0;
}



// Symmetric tridiagonal QL algorithm.
static void tql2(double V[3][3], double d[3], double e[3])
{
    
	const int n=3;
    
	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.
	
	for (int i = 1; i < n; i++) {
		e[i-1] = e[i];
	}
	e[n-1] = 0.0;
	
	double f = 0.0;
	double tst1 = 0.0;
	double eps = pow(2.0,-52.0);
	for (int l = 0; l < n; l++) {
		
		// Find small subdiagonal element
		
		tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
		int m = l;
		while (m < n) {
			if (fabs(e[m]) <= eps*tst1) {
				break;
			}
			m++;
		}
		
		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.
		
		if (m > l) {
			int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)
				
				// Compute implicit shift
				
				double g = d[l];
				double p = (d[l+1] - g) / (2.0 * e[l]);
				double r = hypot2(p,1.0);
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l+1] = e[l] * (p + r);
				double dl1 = d[l+1];
				double h = g - d[l];
				for (int i = l+2; i < n; i++) {
					d[i] -= h;
				}
				f = f + h;
				
				// Implicit QL transformation.
				
				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l+1];
				double s = 0.0;
				double s2 = 0.0;
				for (int i = m-1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot2(p,e[i]);
					e[i+1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i+1] = h + s * (c * g + s * d[i]);
					
					// Accumulate transformation.
					
					for (int k = 0; k < n; k++) {
						h = V[k][i+1];
						V[k][i+1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;
				
				// Check for convergence.
				
			} while (fabs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}
	
	// Sort eigenvalues and corresponding vectors.
	for (int i = 0; i < n-1; i++)
    {
		int k = i;
		double p = d[i];
		for (int j = i+1; j < n; j++)
        {
			if (d[j] < p)
            {
				k = j;
				p = d[j];
			}
		}
		if (k != i)
        {
			d[k] = d[i];
			d[i] = p;
			for (int j = 0; j < n; j++)
            {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
}



const istensor istensor :: random()
{
    itensor t;
    
	const int ndm(3);
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			t(i,j) -= random::uniform(-1.0,1.0);
    
    return symmetricPartOf(t);
}


void istensor :: setRandom()
{
    (*this) = istensor::random();
}



// ivectors holds the three eigenvectors, one per column
void istensor ::  spectralDecomposition(ivector ivectors[3], ivector &evalues) const
{
	const int ndm(3);
	double e[3], d[3], V[3][3];
	for (int i = 0; i < ndm; i++)
		for (int j = 0; j < ndm; j++)
			V[i][j] = (*this)(i,j);
    
	tred2(V, d, e);
	tql2(V, d, e);
	
	for (int i = 0; i < ndm; i++)
    {
		evalues[i] = d[i];
		for (int j = 0; j < ndm; j++)
			ivectors[i][j] = V[j][i];
    }
	
    
	// make sure the eigenvector triad is right-handed.
	// if not, reverse the last eigenvector
	if ( ivector::tripleProduct(ivectors[0], ivectors[1], ivectors[2]) < 0.0)
		for (int i=0; i<ndm; i++) ivectors[2][i] = -ivectors[2][i];
}



const istensor istensor :: deviatoricPart(const istensor& t)
{
    istensor dev(t);
    double   theta3 = (t.trace())/3.0;
    
    dev(0,0) -= theta3;
    dev(1,1) -= theta3;
    dev(2,2) -= theta3;
    
    return dev;
}




const istensor istensor :: symmetricPartOf(const itensor &t)
{
	const int ndm(3);
	istensor s;
	for (int i=0; i<ndm; i++)
	{
		s(i,i) = t(i,i);
		for (int j=i+1; j<ndm; j++)
			s(i,j) = s(j,i) = 0.5*(t(i,j)+t(j,i));
	}
	
	return s;
}



istensor& istensor :: operator=(const istensor &t)
{
	const int ndm(3);
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] = t.a[i][j];
	return *this;
}




istensor& istensor :: operator=(const itensor &t)
{
	*this = symmetricPartOf(t);
	return *this;
}




istensor istensor :: operator-()
{
	const int ndm(3);
	istensor ret;
	
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			ret(i,j) = -a[i][j];
	return ret;
}



istensor&  istensor ::  operator+=(const istensor &t)
{
	const int ndm(3);
    
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] += t.a[i][j];
	return *this;
}


istensor&  istensor ::  operator-=(const istensor &t)
{
	const int ndm(3);
    
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] -= t.a[i][j];
	return *this;
}


istensor&  istensor ::  operator*=(const double alpha)
{
	const int ndm(3);
    
	for (int i=0; i<ndm; i++)
		for (int j=0; j<ndm; j++)
			a[i][j] *= alpha;
	return *this;
}


namespace blue
{
    
    
    istensor operator+(const istensor &t1, const istensor &t2)
    {
        istensor ret(t1);
        ret += t2;
        return ret;
    }
    
    
    itensor   operator+(const istensor &t1, const itensor &t2)
    {
        itensor tt(t1);
        return tt+t2;
    }
    
    
    
    itensor operator+(const itensor &t1, const istensor &t2)
    {
        return t2+t1;
    }
    
    
    
    
    istensor operator-(const istensor &t1, const istensor &t2)
    {
        istensor ret(t1);
        ret -= t2;
        return ret;
    }
    
    
    
    itensor  operator*(const itensor &t1, const istensor &t2)
    {
        itensor ret;
        const int ndm(3);
        
        for (int i=0; i<ndm; i++)
            for (int j=0; j<ndm; j++)
                for (int k=0; k<ndm; k++)
                    ret(i,j) += t1(i,k) * t2.a[k][j];
        
        return ret;
    }
    
    
    itensor  operator*(const istensor &t1, const itensor &t2)
    {
        itensor ret;
        const int ndm(3);
        
        for (int i=0; i<ndm; i++)
            for (int j=0; j<ndm; j++)
                for (int k=0; k<ndm; k++)
                    ret(i,j) += t1(i,k) * t2(k,j);
        
        return ret;
    }
    
    
    
    istensor  operator*(const istensor &t1, const istensor &t2)
    {
        itensor ret;
        const int ndm(3);
        
        for (int i=0; i<ndm; i++)
            for (int j=0; j<ndm; j++)
                for (int k=0; k<ndm; k++)
                    ret(i,j) += t1(i,k) * t2.a[k][j];
        
        return ret;
    }
    
    
    
    ivector  operator*(const istensor &t , const ivector &v)
    {
        ivector w;
        const int ndm(3);
        
        for (int i=0; i<ndm; i++)
            for (int j=0; j<ndm; j++)
                w[i] += t.a[i][j] * v[j];
        return w;
    }
    
    
    istensor  operator*(const istensor &t , const double alpha)
    {
        istensor ret(t);
        ret *= alpha;
        return ret;
    }
    
    
    istensor  operator*(const double alpha, const istensor &t)
    {
        return t*alpha;
    }
    
}


const istensor  istensor :: tensorTimesTensorTransposed(const itensor& F)
{
    istensor FFt;
    
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            FFt(i,j) = 0.0;
            for (int k=0; k<3; k++)
            {
                FFt(i,j) += F(i,k) * F(j,k);
            }
        }
    }
    return FFt;
}



const istensor istensor :: tensorTransposedTimesTensor(const itensor& F)
{
    istensor FtF;
    
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            FtF(i,j) = 0.0;
            for (int k=0; k<3; k++)
            {
                FtF(i,j) += F(k,i) * F(k,j);
            }
        }
    }
    return FtF;
}



const istensor istensor :: FtCF(const itensor& F, const istensor& C)
{
    itensor CF;
    
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            CF(i,j) = 0.0;
            for (int k=0; k<3; k++)
            {
                CF(i,j) += C(i,k) * F(k,j);
            }
        }
    }
    
    istensor mFtCF;
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            mFtCF(i,j) = 0.0;
            for (int k=0; k<3; k++)
            {
                mFtCF(i,j) += F(k,i) * CF(k,j);
            }
        }
    }
    
    return mFtCF;
}



const istensor istensor :: FSFt(const itensor& F, const istensor& S)
{
    itensor SFt;
    
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            SFt(i,j) = 0.0;
            for (int k=0; k<3; k++)
            {
                SFt(i,j) += S(i,k) * F(j,k);
            }
        }
    }
    
    istensor mFSFt;
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            mFSFt(i,j) = 0.0;
            for (int k=0; k<3; k++)
            {
                mFSFt(i,j) += F(i,k) * SFt(k,j);
            }
        }
    }
    
    return mFSFt;
}





irotation :: irotation()
{
	a[0][0] = a[1][1] = a[2][2] = 1.0;
	a[0][1] = a[1][0] = 0.0;
	a[0][2] = a[2][0] = 0.0;
	a[1][2] = a[2][1] = 0.0;
}


irotation :: irotation(const double a, const double b, const double c)
{
	*this = irotation( ivector(a,b,c) );
}


// rotation initialization via exponential map
irotation :: irotation(const ivector& theta)
{
	double norm, f1, f2, f3, sn;
	
	norm = theta.norm();
	f1   = cos(norm);
	
	if (norm < 1.0e-6)
    {
		norm *= norm;
		f2 = 1.0 - norm*(1.0/6.0  - norm*(1.0/120.0 - norm/5040.0));
		f3 = 0.5 - norm*(1.0/24.0 - norm*(1.0/720.0 - norm/40320.0));
    }
	else
    {
		f2    = sin(norm)/norm;
		norm *= 0.5;
		sn    = sin(norm);
		f3    = 0.5 * sn * sn /(norm*norm);
    }
	
	
	a[0][0] =           f1 + f3*theta[0]*theta[0];
	a[0][1] = -theta[2]*f2 + f3*theta[0]*theta[1];
	a[0][2] =  theta[1]*f2 + f3*theta[0]*theta[2];
	
	a[1][0] =  theta[2]*f2 + f3*theta[1]*theta[0];
	a[1][1] =           f1 + f3*theta[1]*theta[1];
	a[1][2] = -theta[0]*f2 + f3*theta[1]*theta[2];
	
	a[2][0] = -theta[1]*f2 + f3*theta[2]*theta[0];
	a[2][1] =  theta[0]*f2 + f3*theta[2]*theta[1];
	a[2][2] =           f1 + f3*theta[2]*theta[2];
}





irotation :: irotation(const iquaternion &quat)
{
	computeFrom(quat);
}


irotation :: irotation(const itensor& t)
{
	for (int i = 0; i<3; i++)
		for (int j=0; j<3; j++)
			this->a[i][j] = t(i,j);
}




void irotation :: beRotationWithoutDrill(const ivector& u)
{
	double f;
	
	if (u[2] > 0.0)
	{
		f   =   1.0 / (1.0 + u[2]);
		a[0][0] =   u[2] + f * u[1] * u[1];
		a[0][1] =        - f * u[1] * u[0];
		a[1][0] =        - f * u[0] * u[1];
		a[1][1] =   u[2] + f * u[0] * u[0];
		a[2][0] = - u[0];
		a[2][1] = - u[1];
	}
    else
	{
		f         =   1.0 / (1.0 - u[2]);
		a[0][0] = - u[2] + f * u[1] * u[1];
		a[0][1] =          f * u[1] * u[0];
		a[1][0] =        - f * u[0] * u[1];
		a[1][1] =   u[2] - f * u[0] * u[0];
		a[2][0] =   u[0];
		a[2][1] = - u[1];
	}
	
	a[0][2] = u[0];
	a[1][2] = u[1];
	a[2][2] = u[2];
}


void irotation :: beRotationWithoutDrill(const ivector &from, const ivector& to)
{
	itensor rot;
	double    dd = from.dot(to);
	ivector   cr = from.cross(to);
	
	rot = dd * identity() + skew( cr ) + 1.0/(1.0+dd)* dyadic(cr,cr);
	*this = rot;
}




void irotation :: computeFrom(const iquaternion &quat)
{
	double q00, q01, q02, q03, q11, q12, q13, q22, q23, q33;
	
	q00 = quat.w()*quat.w()*2.0 - 1.0;
	q01 = quat.w()*quat.x()*2.0;
	q02 = quat.w()*quat.y()*2.0;
	q03 = quat.w()*quat.z()*2.0;
	q11 = quat.x()*quat.x()*2.0;
	q12 = quat.x()*quat.y()*2.0;
	q13 = quat.x()*quat.z()*2.0;
	q22 = quat.y()*quat.y()*2.0;
	q23 = quat.y()*quat.z()*2.0;
	q33 = quat.z()*quat.z()*2.0;
	
	a[0][0] = q00 + q11;
	a[1][0] = q12 + q03;
	a[2][0] = q13 - q02;
	a[0][1] = q12 - q03;
	a[1][1] = q00 + q22;
	a[2][1] = q23 + q01;
	a[0][2] = q13 + q02;
	a[1][2] = q23 - q01;
	a[2][2] = q00 + q33;
}




ivector irotation :: rotationVector() const
{
    ivector v;
    v.extractFrom(*this);
    return v;
}




iquaternion :: iquaternion()
{
	q[0] = q[1] = q[2] = 0.0;
	q[3] = 1.0;
}



iquaternion :: iquaternion(const double q0, const double q1, const double q2, const double q3)
{
	q[0] = q0;
	q[1] = q1;
	q[2] = q2;
	q[3] = q3;
}


iquaternion :: iquaternion(const double* qr)
{
	q[0] = qr[0];
	q[1] = qr[1];
	q[2] = qr[2];
	q[3] = qr[3];
}


iquaternion :: iquaternion(const iquaternion& q_)
{
	q[0] = q_.x();
	q[1] = q_.y();
	q[2] = q_.z();
	q[3] = q_.w();
}


iquaternion :: iquaternion(const irotation &m)
{
	extractFromRotation(m);
}



iquaternion :: iquaternion(const ivector &theta)
{
	// extract the iquaternion such that expmap[theta] = quat->rotationMatrix
	int     i;
    double  thetanr2, rr, fac2;
	
    thetanr2 = theta[0]*theta[0]+theta[1]*theta[1]+theta[2]*theta[2];
    thetanr2 = 0.5*sqrt(thetanr2);
    
    if (thetanr2 < 1e-4)
    {
        rr   = thetanr2*thetanr2;
        fac2 = 0.5 - rr*(840.0 - rr*(42.0- rr))/10080.0;
    }
    else
        fac2 = sin(thetanr2)/thetanr2 * 0.5;
	
    q[3] = cos(thetanr2);
    for (i=0; i<3; i++) q[i] = fac2 * theta[i];
}





// extract iquaternion from a rotation matrix with Spurrier algorithm
void iquaternion :: extractFromRotation(const irotation &mat)
{
	double trace, xm, bigm, dum;
	int    e[]={0,1,2,0,1};
	int    ii, j, k, l, ll, i;
	
	trace = mat.trace();
	ii    = 0;
	xm    = mat(0,0);
	
	if (mat(1,1) > xm)
	{
		xm = mat(1,1);
		ii = 1;
	}
	
	if (mat(2,2) > xm)
	{
		xm = mat(2,2);
		ii  = 2;
	}
	
	bigm = trace > xm ? trace : xm;
	
	if (bigm == trace)
	{
		q[3] = 0.5*sqrt(1.0 + trace);
		dum  = 0.25/ q[3];
		for (i=0; i<3; i++)
		{
			j    = e[i+1];
			k    = e[i+2];
			q[i] = dum*(mat(k,j) - mat(j,k));
		}
	}
	else
	{
		q[ii] = sqrt(0.5*mat(ii,ii) + 0.25*(1.0 - trace));
		dum   = 0.25 / q[ii];
		j     = e[ii+1];
		k     = e[ii+2];
		q[3]  = dum*(mat(k,j) - mat(j,k));
		for (ll=ii+1; ll<=ii+2; ll++)
		{
			l    = e[ll];
			q[l] = dum*(mat(l,ii) + mat(ii,l));
		}
	}
	
	// if, during computations, the iquaternion has lost its orthogonality, normalize it
	double n = norm();
	if ( fabs(n - 1.0) > 1e-6)
	{
		n = 1.0/n;
		for (i=0; i<4; i++) q[i] *= n;
	}
}



iquaternion iquaternion :: conjugate()
{
	return iquaternion(-q[0], -q[1], -q[2], q[3]);
}


double iquaternion :: norm() const
{
	return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
}




void iquaternion :: normalize()
{
	double n(1.0/norm());
	q[0] *= n;
	q[1] *= n;
	q[2] *= n;
	q[3] *= n;
}



void iquaternion :: toIdentity()
{
	q[0] = q[1] = q[2] = 0.0;
	q[3] = 1.0;
}




double& iquaternion :: x()
{
    return q[0];
}

double& iquaternion :: y()
{
    return q[1];
}

double& iquaternion :: z()
{
    return q[2];
}

double& iquaternion :: w()
{
    return q[3];
}


const double& iquaternion :: x() const
{
    return q[0];
}

const double& iquaternion :: y() const
{
    return q[1];
}

const double& iquaternion :: z() const
{
    return q[2];
}

const double& iquaternion :: w() const
{
    return q[3];
}


namespace blue
{
    
    ivector operator*(const iquaternion &quat , const ivector &v)
    {
        irotation m(quat);
        return ivector(m*v);
    }
    
    
    iquaternion& iquaternion :: operator=(const iquaternion &rhs)
    {
        q[0] = rhs.x();
        q[1] = rhs.y();
        q[2] = rhs.z();
        q[3] = rhs.w();
        
        return *this;
    }
    
    
    iquaternion& iquaternion :: operator=(const irotation &rhs)
    {
        *this = iquaternion(rhs);
        return *this;
    }
    
    
    
    iquaternion operator*(const iquaternion &p, const iquaternion &q)
    {
        double pq3(p.x()*q.x()+p.y()*q.y()+p.z()*q.z());
        iquaternion pq(p.w()*q.x() + p.x()*q.w() + p.y()*q.z() - p.z()*q.y(),
                       p.w()*q.y() + p.y()*q.w() + p.z()*q.x() - p.x()*q.z(),
                       p.w()*q.z() + p.z()*q.w() + p.x()*q.y() - p.y()*q.x(),
                       p.w()*q.w() - pq3);
        
        pq.normalize();
        return pq;
    }
    
    std::ostream& operator<<(std::ostream &os, const iquaternion &v)
    {
        os  << std::setw(12) << std::setprecision(8)
        << "[ " << v.x() << " , " << v.y() << " , " << v.z() << " , " << v.w()
        << " ]" << std::flush;
        
        return os;
    }
    
    
    
    
    
    std::ostream& operator<<(std::ostream &os, const vector2 &v)
    {
        os  << std::setw(12) << std::setprecision(8)
        << "[ " << v.x[0] << " , " << v.x[1] << " ]" << std::flush;
        
        return os;
    }
    
}



