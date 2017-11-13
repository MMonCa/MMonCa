/*
 * Author: ignacio.martin@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain
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

#ifndef LKMCLATTICEATOMDIAMONDPARAM_H
#define LKMCLATTICEATOMDIAMONDPARAM_H

#include "LatticeParam.h"
#include "kernel/IDContainer.h"
#include <tcl.h>

namespace LKMC {

class LatticeDiamondParam : public LatticeParam
{
public:
    LatticeDiamondParam(Tcl_Interp *pTcl, const IO::ParameterManager *, const IO::FileParameters *, Kernel::M_TYPE);

    enum PLANE { P100_6=0, P100_7, P100_8, P100_9, P100_10,
    	P110, P111, PNONE };
    enum ELECTROSTATIC_MODEL { GFLS, ATOMISTIC, NONE };
    
    struct SPERParams {
		float _prefSPER[PNONE];
		Kernel::IDContainer<float> _enerSPERFirst;
	};

    SPERParams _sper[Kernel::MAX_IMPURITIES];
    float _twinProb;
	float _shearEffect;

    float _dvpar100_2;
    float _dvpar100_3;
    float _dvpar110;
    float _dvpar111;

    float _dvperp100_2;
    float _dvperp100_3;
    float _dvperp110;
    float _dvperp111;

    // charged configurations
    double              _E_M0;
    double              _E_P0;
    double              _g_M;
    double              _g_P;
    double              _g_0;
    ELECTROSTATIC_MODEL _electrostaticModel;

    //epitaxy
    bool  _model_simplified;

    //the other parameters are fixed numbers and depend only on the number of them.
    struct EpiParams
    {
    	Kernel::IDContainer<float> _formationFirst;  //formation first depends on the cluster configuration (Si2Ge1, etc...)
        float _formationSecond[Kernel::MAX_IMPURITIES];
    	float _formationThird [Kernel::MAX_IMPURITIES];
    	Kernel::P_TYPE _pairPrecursor;  //when adsorption by pairs is faster

		float _migrationF;
		float _migrationS;
		float _barrierEpi;
		float _barrierPairEpi;
		float _barrierPrecursor;
		float _barrierDes;
		float _prefEpi;
		float _prefEtch;
		float _prefMig;
		float _prefDes;
		float _speedUpRatio;
    };
    EpiParams _epi[Kernel::MAX_IMPURITIES];

private:
    void readSPER   (Tcl_Interp *pTcl, const IO::ParameterManager *, const IO::FileParameters *, Kernel::M_TYPE);
    void readEpitaxy(Tcl_Interp *pTcl, const IO::ParameterManager *, const IO::FileParameters *, Kernel::M_TYPE);
};

}

#endif
