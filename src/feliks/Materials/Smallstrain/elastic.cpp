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



#include "elastic.h"

#include <cmath>
#include "Math/statistics.h"

#include "Io/usercommand.h"
#include "Io/logger.h"
#include "Io/message.h"
#include "Math/tensor.h"

using namespace blue;


elasticMaterial :: elasticMaterial(const std::string& name,
			   const double _reft,
               const double _thermal,
               const double _E,
               const double _nu,
               const std::map<string, std::pair<double, double> > & _eigenstrains)
:
    smallStrainMaterial(_eigenstrains),
    E(_E), nu(_nu), thermal(_thermal), temp0(_reft)
{
    lambda = E*nu/(1.0-2.0*nu)/(1.0+nu);
    mu     = E/2.0/(1.0+nu);
    bulk   = lambda + 2.0/3.0 * mu;

  /*  if (rho > 0.0)
    {
    	 cp     = sqrt((lambda+2.0*mu)/rho);
    	 cs     = sqrt(2.0*mu/rho);
    }*/

}




bool elasticMaterial :: check() const
{
	if (mu > 0 && lambda+2.0*mu > 0) return true;
	else return false;
}




/* an object of the type "elastic" creates a material point of type elastic. The
 material point holds information that is not part of the material itself but
 that is particular of the specific (physical) point.
 */
smallStrainMP* elasticMaterial :: createMaterialPoint(const std::vector<double>* comp) const
{
	smallStrainMP* mp = new elasticMP(*this);
	return mp;
}




double elasticMaterial :: density() const
{
    return rho;
}






/* this function is always called once the material is defined, so apart from
 printing its information, we take the opportunity to clean up some of its
 data, in particular, setting all the possible constants
 */
void elasticMaterial :: print(std::ostream &of) const
{
    of  << "\n   Small strain, elastic, isotropic material ";
    of  << "\n   Young modulus   E      : " << E;
    of  << "\n   Poisson ratio   nu     : " << nu;
    of  << "\n   Lame constants  lambda : " << lambda;
    of  << "\n                   mu     : " << mu;
    of  << "\n   Bulk modulus    k      : " << bulk;
    of  << "\n   Density                : " << rho;
    of  << "\n   Reference temperature  : " << temp0;
    if (rho > 0.0)
    {
        of  << "\n   Wave velocities C_p    : " << cp;
        of  << "\n                   C_s    : " << cs;
    }
    of << std::endl;
}




void elasticMaterial :: setRandom()
{
    E      = random::uniform(1.0, 10000.0);
    nu     = random::uniform(0.05, 0.45);
    rho    = random::uniform(1.0, 100.0);
    lambda = E*nu/(1.0-2.0*nu)/(1.0+nu);
    mu     = E/2.0/(1.0+nu);
    cp     = sqrt((lambda+2.0*mu)/rho);
    cs     = sqrt(2.0*mu/rho);
    bulk   = lambda + 2.0/3.0 * mu;
}




bool elasticMaterial :: test()
{
    bool isok = true;
    setRandom();
    smallStrainMP* p = this->createMaterialPoint();
    
    DebugMessage("");
    DebugMessage("Testing small strain elastic material");
    isok = p->testImplementation();
    if (!isok) std::cout << "\n Error in test" << std::endl;
    
    return isok;
}




elasticMP :: elasticMP(const elasticMaterial &m) :
smallStrainMP(m),
theElasticMaterial(m),
tn(m.temp0), tc(m.temp0)
{
    En.setZero();
    Ec.setZero();
}




elasticMP :: ~elasticMP()
{
    
}




void elasticMP :: commitCurrentState()
{
    tn    = tc;
    En    = Ec;
    tempn = tempc;
}




/* Given the fourth order tensor of elasticities C, and two vectors v, w
 * compute the second order tensor T with components
 *   T_ij = C_ipjq v_p w_q
 *
 *  which, for linear isotropic materials is just
 *  T = lambda gradV otimes gradW + mu gradW otimes gradV + mu (gradV * gradW) One
 *
 *  Note that the result is not symmetric.
 */

void elasticMP :: contractWithTangent(const ivector &v1, const ivector &v2, itensor &T) const
{
	itensor A, B;
	
	T = itensor::identity();
	T *= v1.dot(v2) * theElasticMaterial.mu;
	
	A = itensor::dyadic(v1, v2);
	B = A.transpose();
	
	A *= theElasticMaterial.lambda;
	B *= theElasticMaterial.mu;
	
	T += A;
	T += B;
}




void elasticMP :: contractWithDeviatoricTangent(const ivector &v1, const ivector &v2, itensor &T) const
{
	const double&  mu = theElasticMaterial.mu;
    
	T = mu * v1.dot(v2) * itensor::identity()
    +   mu *              itensor::dyadic(v2, v1)
    -   2.0/3.0*mu *      itensor::dyadic(v1, v2);
}




double elasticMP :: energyDissipationPotential() const
{
    return storedEnergy();
}




void elasticMP :: materialTangentTimesSymmetricTensor(const istensor& M, istensor& CM) const
{
    CM = theElasticMaterial.lambda * M.trace() * istensor::identity()
    +    2.0* theElasticMaterial.mu * M;
}




double elasticMP :: plasticSlip() const
{
    return 0.0;
}



void elasticMP :: resetCurrentState()
{
    tc    = tn;
    Ec    = En;
    tempc = tempn;
}



double elasticMP :: storedEnergy() const
{
    const double  tr = Ec.trace();
    return  theElasticMaterial.mu*Ec.squaredNorm() + 0.5*theElasticMaterial.lambda*tr*tr;
}




void elasticMP :: stress(istensor& sigma) const
{
	sigma = (2.0*theElasticMaterial.mu)*Ec +
        (theElasticMaterial.lambda*Ec.trace() -
        		3.0*theElasticMaterial.bulk*theElasticMaterial.thermal*(tempc-theElasticMaterial.temp0)) * istensor::identity();
}




void elasticMP :: updateCurrentState(const double theTime, const istensor& strain, const double temp)
{
    tc    = theTime;
    Ec    = strain;
    tempc = temp;
}




double elasticMP :: volumetricStiffness() const
{
	return theElasticMaterial.bulk;
}

