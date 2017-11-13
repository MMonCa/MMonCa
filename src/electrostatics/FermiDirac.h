/*
 * FermiDirac.h
 *
 *  Created on: Nov 21, 2014
 *      Author: ignacio.martin@imdea.org
 */

#ifndef FERMIDIRAC_H_
#define FERMIDIRAC_H_

namespace Electrostatics {

class FermiDirac {
public:
	static double phalf(double x);
	static double mhalf(double x);
};

} /* namespace Electrostatics */
#endif /* FERMIDIRAC_H_ */
