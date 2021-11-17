/*
 * FermiDirac.cpp
 *
 *  Created on: Nov 21, 2014
 *      Author: ignacio.martin@imdea.org
 */

#include "FermiDirac.h"
#include <cmath>

namespace Electrostatics {

//Based in the model published in D. Bednarczyk and J. Bednarczyk, Phys. Lett. A, 64, 409 (1978)
//but the expressions are from equations 22 to 24 of J. S. Blakemore, Solid-State Electronics, 25, 1067 (1982)
double FermiDirac::phalf(double x)
{
	if (x<-700)
		return 0;
	 double mu = x*x*x*x + 50 +
			 x* 33.6 * ( 1. - .68 * exp( -0.17 * (x + 1)*(x+1)));
	 double xn = 3. * sqrt(M_PI) / (4* pow(mu, 3./8.));
	 return 1./(exp(- x ) + xn );
}

//Based in equations 6 to 7 of X. Aymerich-Humet, F. Serra-Mestres, and J. Millan, J. Appl. Phys., 54, 2850 (1983)
double FermiDirac::mhalf(double x)
{
	if (x<-700)
		return 0;
	 const double a =  sqrt( 1. + 15. / 8. + 1 / 160.);
	 const double b = 1.8 + 0.61 / 2.;
	 const double c = 2. + ( 2 - sqrt( 2.) ) * sqrt(2.);
	 return 1./( 1 / sqrt( b + x + pow(( pow(fabs( x - b ), c) + pow(a,c)), ( 1./c))  ) +
			 exp(-x ) / sqrt(M_PI) ) / sqrt(M_PI);
}

} /* namespace Electrostatics */
