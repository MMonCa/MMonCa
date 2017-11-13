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

#include "LatticeDiamondParam.h"
#include "EpiGasParam.h"
#include "io/ParameterManager.h"
#include "io/FileParameters.h"

using std::string;
using std::map;
using std::vector;

namespace LKMC {

LatticeDiamondParam::LatticeDiamondParam(Tcl_Interp *pTcl, const IO::ParameterManager *pPM, const IO::FileParameters *pPar, Kernel::M_TYPE mt)
: LatticeParam(pPM, pPar, mt)
{
	_type = (pPM->getMaterial(mt)._binary? DIAMOND2 : DIAMOND);

	readSPER   (pTcl, pPM, pPar, mt);
	readEpitaxy(pTcl, pPM, pPar, mt);
}

void LatticeDiamondParam::readSPER(Tcl_Interp *pTcl, const IO::ParameterManager *pPM, const IO::FileParameters *pPar, Kernel::M_TYPE mt)
{
	string base       = pPM->getMaterialName(mt) + "/SPER/";

	//SPER
	_shearEffect      = pPar->getFloat(base + "shear.effect");
	_twinProb          = pPar->getFloat(base + "twin.probability");

    _dvpar100_2        = pPar->getFloat(base + "dvpar.100.2");
    _dvpar100_3        = pPar->getFloat(base + "dvpar.100.3");
    _dvpar110          = pPar->getFloat(base + "dvpar.110");
    _dvpar111          = pPar->getFloat(base + "dvpar.111");

    _dvperp100_2       = pPar->getFloat(base + "dvperp.100.2");
    _dvperp100_3       = pPar->getFloat(base + "dvperp.100.3");
    _dvperp110         = pPar->getFloat(base + "dvperp.110");
    _dvperp111         = pPar->getFloat(base + "dvperp.111");

    _E_M0              = pPar->getFloat(base + "energy.M0");
    _E_P0              = pPar->getFloat(base + "energy.P0");
    _g_M               = pPar->getFloat(base + "degeneracy.M");
    _g_P               = pPar->getFloat(base + "degeneracy.P");
    _g_0               = pPar->getFloat(base + "degeneracy.0");

    if (pPar->getString(base + "electrostatic.model") == "GFLS")
    	_electrostaticModel = GFLS;
    else if (pPar->getString(base + "electrostatic.model") == "atomistic")
    	_electrostaticModel = ATOMISTIC;
    else
    	_electrostaticModel = NONE;

    map<string,float> energy1Map = pPar->getFloatMap(base + "energy.1st");

    for(map<string,float>::iterator it=energy1Map.begin(); it != energy1Map.end(); ++it)
	{
		vector<string> tokens;
		IO::ParameterManager::getTokens(it->first, '+', tokens);
		if(tokens.size() != 2)
			ERRORMSG(base + "energy.1st/" << it->first << " syntax error. Two tokens and the '+' sign are required.");

		Kernel::P_TYPE pt = pPM->getParticleNumber(mt, tokens[0]);
		pt = pPM->getFamily(pt);
		if(pt == Kernel::UNDEFINED_TYPE)
			ERRORMSG("Wrong particle index " << base+"energy.1st" << ": " << it->first);

		Kernel::ID clusterID = pPM->getEpiID(mt, tokens[1]);
		if(clusterID._pt.empty())
			ERRORMSG("Wrong ID in " << base+"energy.1st" << ": " << it->first);
		clusterID._pos.clear();
		if(tokens[1] != pPM->getEpiIDName(clusterID))
			WARNINGMSG(base + "energy.1st" << ": Suggested syntax for '" << tokens[1] << "' is '" << pPM->getEpiIDName(clusterID) << "'");

		unsigned hash = _sper[pt]._enerSPERFirst.addHash(clusterID);
		//Now I have the list of clusters, so I can ask for their potentials...
		_sper[pt]._enerSPERFirst._map[hash]= it->second;
	}

    //read string - float maps
   string names[7] = { "prefactor.100.6", "prefactor.100.7", "prefactor.100.8", "prefactor.100.9",
		   "prefactor.100.10", "prefactor.110", "prefactor.111" };
   for(unsigned i=0; i < 7; ++i)
   {
		map<string, float> theMap = pPar->getFloatMap(base+names[i]);
		for(map<string, float>::iterator it=theMap.begin(); it!=theMap.end(); ++it)
		{
			Kernel::P_TYPE pt = pPM->getParticleNumber(mt, it->first);
			if(pt == Kernel::UNDEFINED_TYPE)
				ERRORMSG("Wrong particle index " << base+names[i] << ": " << it->first);
			pt = pPM->getFamily(pt);
			switch(i)
			{
			case 0:
				_sper[pt]._prefSPER[P100_6] = it->second; break;
			case 1:
				_sper[pt]._prefSPER[P100_7] = it->second; break;
			case 2:
				_sper[pt]._prefSPER[P100_8] = it->second; break;
			case 3:
				_sper[pt]._prefSPER[P100_9] = it->second; break;
			case 4:
				_sper[pt]._prefSPER[P100_10] = it->second; break;
			case 5:
				_sper[pt]._prefSPER[P110] = it->second; break;
			case 6:
				_sper[pt]._prefSPER[P111] = it->second; break;
			default:
				ERRORMSG("Wrong case in LatticeDiamondParam line " << __LINE__);
			}
		}
   }
}

void LatticeDiamondParam::readEpitaxy(Tcl_Interp *pTcl, const IO::ParameterManager *pPM, const IO::FileParameters *pPar, Kernel::M_TYPE mt)
{
    string base       = Domains::global()->PM()->getMaterialName(mt) + "/Epitaxy";

    _model_simplified     = pPar->getBool (base+"/model.simplified");

    map<string,float> frm1Map = pPar->getFloatMap(base + "/formation.1st");
    map<string,float> frm2Map = pPar->getFloatMap(base + "/formation.2nd");
    map<string,float> frm3Map = pPar->getFloatMap(base + "/formation.3rd");
    //read the flexible introduction of precursor prefactors depending on pressure and temperature
	pPar->loadProcedure(pTcl, base+"/prefactor.precursor", MAX_EPI_GASES*2 + 2);

    for(map<string,float>::iterator it=frm1Map.begin(); it != frm1Map.end(); ++it)
    {
    	vector<string> tokens;
    	IO::ParameterManager::getTokens(it->first, '+', tokens);
    	if(tokens.size() != 2)
    		ERRORMSG(base + "/formation.1st/" << it->first << " syntax error. Two tokens and the '+' sign are required.");

    	Kernel::P_TYPE pt = pPM->getParticleNumber(mt, tokens[0]);
		pt = pPM->getFamily(pt);
		if(pt == Kernel::UNDEFINED_TYPE)
			ERRORMSG("Wrong particle index " << base+"/formation.1st" << ": " << it->first);

		Kernel::ID clusterID = pPM->getEpiID(mt, tokens[1]);
		if(clusterID._pt.empty())
			ERRORMSG("Wrong ID in " << base+"/formation.1st" << ": " << it->first);
		clusterID._pos.clear();
		if(tokens[1] != pPM->getEpiIDName(clusterID))
			WARNINGMSG(base + "/formation.1st" << ": Suggested syntax for '" << tokens[1] << "' is '" << pPM->getEpiIDName(clusterID) << "'");

		unsigned hash = _epi[pt]._formationFirst.addHash(clusterID);
		//Now I have the list of clusters, so I can ask for their potentials...
		_epi[pt]._formationFirst._map[hash]= it->second;
    }

    for(Kernel::P_TYPE pt1 = 0; pt1 < Kernel::MAX_IMPURITIES; ++pt1)
    {
    	_epi[pt1]._barrierEpi = _epi[pt1]._migrationF = _epi[pt1]._migrationS = _epi[pt1]._barrierPrecursor = _epi[pt1]._barrierDes = 5;
    	_epi[pt1]._prefEpi    = _epi[pt1]._prefEtch   = _epi[pt1]._prefMig    = _epi[pt1]._prefDes = 0;
    	for(Kernel::P_TYPE pt2 = 0; pt2 < Kernel::MAX_IMPURITIES; ++pt2)
    		_epi[pt1]._formationSecond[pt2] = _epi[pt1]._formationThird[pt2] = 0; //init values
    }

    for(map<string,float>::iterator it=frm2Map.begin(); it != frm2Map.end(); ++it)
	{
    	vector<string> tokens;
    	IO::ParameterManager::getTokens(it->first, '+', tokens);
    	if(tokens.size() != 2)
    		ERRORMSG(base + "/formation.2nd/" << it->first << " syntax error. Two tokens and the '+' sign are required.");

    	Kernel::P_TYPE pt1 = pPM->getParticleNumber(mt, tokens[0]);
		Kernel::P_TYPE pt2 = pPM->getParticleNumber(mt, tokens[1]);
		if(pt1 == Kernel::UNDEFINED_TYPE || pt2 == Kernel::UNDEFINED_TYPE)
			ERRORMSG("Wrong particle index " << base+"/formation.2nd" << ": " << it->first);

		pt1 = pPM->getFamily(pt1);
		pt2 = pPM->getFamily(pt2);
		_epi[pt1]._formationSecond[pt2] = it->second;
	}

    for(map<string,float>::iterator it=frm3Map.begin(); it != frm3Map.end(); ++it)
	{
    	vector<string> tokens;
		IO::ParameterManager::getTokens(it->first, '+', tokens);
		if(tokens.size() != 2)
			ERRORMSG(base + "/formation.3rd/" << it->first << " syntax error. Two tokens and the '+' sign are required.");

		Kernel::P_TYPE pt1 = pPM->getParticleNumber(mt, tokens[0]);
		Kernel::P_TYPE pt2 = pPM->getParticleNumber(mt, tokens[1]);
		if(pt1 == Kernel::UNDEFINED_TYPE || pt2 == Kernel::UNDEFINED_TYPE)
			ERRORMSG("Wrong particle index " << base+"/formation.3st" << ": " << it->first);

		pt1 = pPM->getFamily(pt1);
		pt2 = pPM->getFamily(pt2);
		_epi[pt1]._formationThird[pt2] = it->second;
	}

    //read string - float maps
    string names[11] = { "/migration.1st", "/migration.2nd", "/barrier.epi", "/prefactor.epi", "/prefactor.etch", "/prefactor.mig",
    		"/barrier.precursor", "/prefactor.desorption", "/barrier.desorption", "/speedup.factor", "/barrier.pair" };
    for(unsigned i=0; i < 11; ++i)
    {
    	map<string, float> theMap = pPar->getFloatMap(base+names[i]);
		for(map<string, float>::iterator it=theMap.begin(); it!=theMap.end(); ++it)
		{
			Kernel::P_TYPE pt = pPM->getParticleNumber(mt, it->first);
			if(pt == Kernel::UNDEFINED_TYPE)
				ERRORMSG("Wrong particle index " << base+names[i] << ": " << it->first);
			pt = pPM->getFamily(pt);
			switch(i)
			{
			case 0:
				_epi[pt]._migrationF = it->second; break;
			case 1:
				_epi[pt]._migrationS = it->second; break;
			case 2:
				_epi[pt]._barrierEpi = it->second; break;
			case 3:
				_epi[pt]._prefEpi = it->second; break;
			case 4:
				_epi[pt]._prefEtch = it->second; break;
			case 5:
				_epi[pt]._prefMig = it->second; break;
			case 6:
				_epi[pt]._barrierPrecursor = it->second; break;
			case 7:
				_epi[pt]._prefDes = it->second; break;
			case 8:
				_epi[pt]._barrierDes = it->second; break;
			case 9:
				_epi[pt]._speedUpRatio = it->second; break;
			case 10:
				_epi[pt]._barrierPairEpi = it->second; break;
			default:
				ERRORMSG("Wrong case in LatticeDiamondParam line " << __LINE__);
			}
		}
    }

    //read string, string maps
    for(Kernel::P_TYPE pt = 0; pt < Kernel::MAX_IMPURITIES; ++pt)
    	_epi[pt]._pairPrecursor = Kernel::UNDEFINED_TYPE; //init
    map<string, string> strstrMap = pPar->getStringMap(base+"/pair.epi");
	for(map<string, string>::iterator it=strstrMap.begin(); it!=strstrMap.end(); ++it)
	{
		Kernel::P_TYPE pt[2] = { pPM->getParticleNumber(mt, it->first) , pPM->getParticleNumber(mt, it->second) };
		if(pt[0] == Kernel::UNDEFINED_TYPE || pt[1] == Kernel::UNDEFINED_TYPE)
			ERRORMSG("Wrong particle index " << base+"/pair.epi" << ": " << it->first << " or " << it->second);
		pt[0] = pPM->getFamily(pt[0]);
		pt[1] = pPM->getFamily(pt[1]);
		_epi[pt[0]]._pairPrecursor = pt[1];
	}
}

}

