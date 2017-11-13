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
/* description of elastic materials
 *
 * i. romero, august 2002
 */

#ifndef _elastic_h
#define _elastic_h


#include "smallstrain.h"
#include "Materials/material.h"
#include "Math/tensor.h"

class commandLine;
class elasticMP;



class elasticMaterial : public smallStrainMaterial
{
public:
                            elasticMaterial(const std::string& name,
                            		        const double _reft,
                                            const double _thermal,
                                            const double _E,
                                            const double _nu,
                                            const std::map<std::string, std::pair<double, double> > & _eigenstrains);
	
	virtual bool            check() const;
    virtual smallStrainMP*  createMaterialPoint(const std::vector<double>* comp=0) const;
    virtual double          density() const;
	virtual void            print(std::ostream &of=std::cout) const;
    virtual void            setRandom();
    virtual bool            test();

    
private:	
    double                  E, nu, thermal, temp0;
	double                  lambda, mu;
	double                  bulk;
	double                  cp, cs;
    double                  rho;

	friend class            elasticMP;
    friend class            alloyMP;
};










class elasticMP : public smallStrainMP
{
public:
                            elasticMP(const elasticMaterial &m);
    virtual                 ~elasticMP();
    
    
    // three dimensional response
    virtual double          storedEnergy() const;
    virtual double          energyDissipationPotential() const;
    virtual void            stress(blue::istensor& sigma) const;
    virtual double          plasticSlip() const;
    

	virtual void            contractWithTangent(const blue::ivector &v1, const blue::ivector &v2, blue::itensor &T) const;
	virtual void            contractWithDeviatoricTangent(const blue::ivector &v1, const blue::ivector &v2, blue::itensor &T) const;
    virtual void            materialTangentTimesSymmetricTensor(const blue::istensor& M, blue::istensor& CM) const;
	virtual double          volumetricStiffness() const;
    
        
    // bookkeeping
    virtual void            updateCurrentState(const double theTime, const blue::istensor& strain, const double temp);
    virtual void            commitCurrentState();
    virtual void            resetCurrentState();

    virtual void            setRandom(){}


private:
	const elasticMaterial&   theElasticMaterial;
    double                  tn, tc;
    blue::istensor                En, Ec;
    double                  tempn, tempc;
};



#endif
