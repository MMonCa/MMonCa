/*
 * ExtractCmd.cpp
 *
 *  Created on: Feb 15, 2011
 *
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

#include "ExtractCmd.h"
#include "FileParameters.h"
#include "domains/SimData.h"
#include "domains/MCClient.h"
#include "kernel/Domain.h"
#include "kernel/MeshElement.h"
#include "kernel/Coordinates.h"
#include "io/ParameterManager.h"
#include "kernel/Domain.h"
#include "okmc/Defect.h"
#include <sstream>
#include <fstream>

namespace IO
{

ExtractCmd::ExtractCmd(Tcl_Interp *p, int argc, const char *argv[]) : Command(p, argc, argv)
{
}

int ExtractCmd::operator()()
{
	std::stringstream ss;

	if(specified("dimension"))
	{
		unsigned dim = getInt("dimension");
		Domains::global()->data()->setDimension(dim);
	}

	Kernel::Coordinates m;
	Kernel::Coordinates M;

	Domains::global()->client()->getCell(m, M);

	if (specified("min.x"))	m._x = getFloat("min.x");
	if (specified("min.y"))	m._y = getFloat("min.y");
	if (specified("min.z"))	m._z = getFloat("min.z");
	if (specified("max.x"))	M._x = getFloat("max.x");
	if (specified("max.y"))	M._y = getFloat("max.y");
	if (specified("max.z"))	M._z = getFloat("max.z");

	Domains::global()->data()->setMinMax(m,M);  //for collapse only.

	if(specified("ac.mean"))
	{
		Kernel::Coordinates mean;

		Domains::global()->data()->getACInterfaceMean(m, M, mean);
		ss << mean._x << " " << mean._y << " " << mean._z;
	}
	else if(specified("ac.stdev"))
	{
		Kernel::Coordinates stdev;

		Domains::global()->data()->getACInterfaceStdev(m, M, stdev);
		ss << stdev._x << " " << stdev._y << " " << stdev._z;
	}
	else if(specified("ac.min"))
	{
		Kernel::Coordinates acmin, acmax;

		Domains::global()->data()->getACInterfaceMinMax(m, M, acmin, acmax);
		ss << acmin._x << " " << acmin._y << " " << acmin._z;
	}
	else if(specified("ac.max"))
	{
		Kernel::Coordinates acmin, acmax;

		Domains::global()->data()->getACInterfaceMinMax(m, M, acmin, acmax);
		ss << acmax._x << " " << acmax._y << " " << acmax._z;
	}
	else if(specified("ac.coverage"))
	{
		ss << Domains::global()->data()->getACInterfaceCoverage(m, M);
	}
	else if(specified("configuration"))
	{
		std::string filename;
		if(specified("filename")) {
			filename = getString("filename");
		}
		Domains::global()->getFileParameters()->dump(ss, filename);
	}
	else if(specified("diffusivity"))
	{
		std::string name = getString("name");
		std::string mate = getString("material");
		bool bMacro = specified("macroscopic");

		ss << Domains::global()->data()->getDiffusivity(mate, name, bMacro);
	}
	else if(specified("jumps"))
	{
		std::string name = getString("name");
		ss << Domains::global()->data()->getJumps(name);
	}
	else if(specified("profile.mobile"))
		// transforms jumps into concentration by using [X] = #jumps / (Delta_time(s) * Volume(cm^3) * jump_freq(s^-1) )
	{
		std::string mate = specified("material")?   getString("material") : "";
		Kernel::M_TYPE mt = Kernel::MAX_MATERIALS;
		if(mate.size())
		{
			mt = Domains::global()->PM()->getMaterialNumber(mate);
			if(mt == Kernel::MAX_MATERIALS)
				ERRORMSG("profile.mobile: Material name '" << mate << "'not recognized");
		}
		if(specified("state"))
		{
			std::string name = getString("name");
			std::string st = getString("state");
			ss << Domains::global()->data()->getProfileMobile(name,st, mate);
		}
		else
		{
			std::string name = getString("name");
			ss << Domains::global()->data()->getProfileMobile(name,"", mate);
		}

	}
	else if(specified("reset"))
		Domains::global()->data()->reset();
	else if(specified("profile"))
		//allows alloys
	{
		std::string mate = specified("material") ? getString("material") : "";
		std::string name;
		if(specified("name"))
			name = getString("name");
		else if(specified("material"))
			name = "";                 // This extract the profile of the A atoms of the AB alloy
		else
			ERRORMSG("profile: At least 'material' or 'name' must be specifed");
		Kernel::M_TYPE mt = Kernel::MAX_MATERIALS;
		std::string st;
		if(mate.size())
		{
			mt = Domains::global()->PM()->getMaterialNumber(mate);
			if(mt == Kernel::MAX_MATERIALS)
				ERRORMSG("profile: Material name '" << mate << "'not recognized");
		}
		if (name == "lkmc.defects" || name == "lkmc.ac") {
			ss << Domains::global()->data()->getLKMCProfile(name);
		}
		else if (name == "A.atoms" || name == "B.atoms")
		{
			ss << Domains::global()->data()->getAtomProfile(name, mate);
		}
		else
		{
			std::string defe = specified("defect") ? getString("defect") : "";
			std::string st = specified("state") ? getString("state") : "";
			ss << Domains::global()->data()->getParticleProfile(name, defe, st, mate);
		}
	}
	else if(specified("defects"))
	{
		std::vector<Domains::ParticleData> data;
		Domains::global()->data()->getOKMCParticles(data, getStrings("defects", 1));
		                    //      mt        coord              part    def  ID           defectID
//		typedef boost::tuple<Kernel::M_TYPE, Kernel::Coordinates, int,   int, std::string, unsigned> ParticleData;
		unsigned nParticles = Domains::global()->PM()->getNParticles();
		for(std::vector<Domains::ParticleData>::iterator it=data.begin(); it!=data.end(); ++it)
			ss << Domains::global()->PM()->getParticleName(it->get<0>(), it->get<2>() % nParticles ) << " " <<
			it->get<1>() << " " << OKMC::Defect::getEName(Kernel::Event::E_TYPE(it->get<3>())) << " " << it->get<4>() << " " <<
			it->get<5>() << std::endl;
	}
	else if(specified("strain"))
	{
		std::string name = getString("name");
		ss << Domains::global()->data()->getStrain(name);
	}
	else if(specified("stress"))
	{
		std::string name = getString("name");
		ss << Domains::global()->data()->getStress(name);
	}
	else if(specified("electrostatic"))
	{
		std::string name = getString("name");
		ss << Domains::global()->data()->getElectrostatic(name);
	}
	else if(specified("time"))
		ss << Domains::global()->getTime();
	else if(specified("count.particles") || specified("count.defects"))
		//allows alloys
	{
		std::string part = specified("particle")?   getString("particle") : "";
		std::string defe = specified("defect")?     getString("defect")   : "";
		std::string sID  = specified("ID")?         getString("ID")       : "";
		std::string mate = specified("material")?   getString("material") : "";
		std::string stat = specified("state")?		getString("state")    : "";
		unsigned minSize = specified("min.size")?   getInt("min.size") : 0;
		Kernel::M_TYPE mt = Kernel::MAX_MATERIALS;
		if(mate.size())
		{
			mt = Domains::global()->PM()->getMaterialNumber(mate);
			if(mt == Kernel::MAX_MATERIALS)
				ERRORMSG("count: Material name '" << mate << "'not recognized");
		}
		if(specified("count.defects") && !stat.empty())
			ERRORMSG("Sorry, states in defects not implemented yet. Please, try with particles.");

		if(specified("count.particles"))
			ss << Domains::global()->data()->count(mt, part, defe, sID, minSize, false, stat); //false = particles
		else
			ss << Domains::global()->data()->count(mt, part, defe, sID, minSize, true, ""); //true = defects

	}
	else if(specified("count.positions"))
	{
		std::string mate = specified("material")?   getString("material") : "";
		std::string pos  = getString("position");
		std::string defe = specified("defect")?     getString("defect")   : "";
		std::string sID  = specified("ID")?         getString("ID")       : "";
		Kernel::M_TYPE mt = Kernel::MAX_MATERIALS;
		if(mate.size())
		{
			mt = Domains::global()->PM()->getMaterialNumber(mate);
			if(mt == Kernel::MAX_MATERIALS)
				ERRORMSG("count: Material name '" << mate << "'not recognized");
		}
		ss << Domains::global()->data()->countPos(mt, pos, defe, sID);
	}
	else if(specified("defect.radius"))
	{
		std::string defe = specified("defect")?     getString("defect")   : "";
		std::string sID  = specified("ID")?         getString("ID")       : "";
		std::string mate = specified("material")?   getString("material") : "";
		unsigned minRad = specified("min.radius")?  getFloat("min.radius") : 0;
		Kernel::M_TYPE mt = Kernel::MAX_MATERIALS;
		if(mate.size())
		{
			mt = Domains::global()->PM()->getMaterialNumber(mate);
			if(mt == Kernel::MAX_MATERIALS)
				ERRORMSG("defect.radius: Material name '" << mate << "'not recognized");
		}
		ss << Domains::global()->data()->radius(mt, defe, sID, minRad);
	}
	else if(specified("histogram"))
	{
		std::string defe = getString("defect");
		std::string mate = getString("material");
		Kernel::M_TYPE mt = Domains::global()->PM()->getMaterialNumber(mate);
		if(mt == Kernel::MAX_MATERIALS)
			ERRORMSG("histogram: Material name '" << mate << "'not recognized");
		ss << Domains::global()->data()->getHistogram(mt, defe);
	}
	else if(specified("coordination"))
	{
		double r = getFloat("radius");
		Kernel::Coordinates c = getCoordinates("coord");
		int typ = specified("type") ? getInt("type") : -1;
		ss << Domains::global()->data()->getCoordination(r, c, typ);
	}
	else if(specified("mixing.enthalpy"))
	{
		std::string mate = specified("material")?   getString("material") : "";
		Kernel::M_TYPE mt = Kernel::MAX_MATERIALS;
		std::vector<float> ov;
		if(mate.size())
		{
			mt = Domains::global()->PM()->getMaterialNumber(mate);
			if(mt == Kernel::MAX_MATERIALS)
				ERRORMSG("mixing.enthalpy: Material name '" << mate << "'not recognized");
		}
		ov = Domains::global()->data()->getMixEnthalpy(mate);
		float x = 0;
		for( std::vector<float>::iterator it = ov.begin(); it != ov.end(); ++it)
		{
			ss << x << " " << *it << std::endl;
			x += 0.001;
		}
	}
	else if(specified("count.displaced"))
	{
		std::string mate = specified("material")?   getString("material") : "";
		if(mate.empty())
			ERRORMSG("count.displaced needs the material to be specified.");
		Kernel::M_TYPE mt = Kernel::MAX_MATERIALS;
		mt = Domains::global()->PM()->getMaterialNumber(mate);
		if(mt == Kernel::MAX_MATERIALS)
			ERRORMSG("amorphous.fraction: Material name '" << mate << "'not recognized");
		const bool bAmorph = specified("amorphous.active") ? getBool("amorphous.active") : true;
		const float threshold = specified("threshold") ? getFloat("threshold") : -1;
		ss << Domains::global()->data()->countDisplaced(mate, bAmorph,threshold);
	}
	else if(specified("amorphous.fraction"))
	{
		ss << Domains::global()->data()->amorphousFraction();
	}
	else if(specified("profile.damage"))
	{
		ss << Domains::global()->data()->getDamageProfile();
	}
	else if(specified("fuzz"))
	{
		ss << Domains::global()->data()->getFuzz();
	}    
        else if(specified("dose"))
	{		
		std::string mate = getString("material");
		Kernel::M_TYPE mt = Domains::global()->PM()->getMaterialNumber(mate);
		if(mt == Kernel::MAX_MATERIALS)
			ERRORMSG("dose: Material name '" << mate << "'not recognized");
		ss << Domains::global()->data()->getDose(mt);
	}
        else if(specified("balance"))
        {
                std::string mate = getString("material");
                std::string name = getString("name");	                
                ss << Domains::global()->data()->getBalance(name, mate);
        }
	else
		ERRORMSG("No known option specified");

	if(specified("file"))
	{
		std::string filename = getString("file");
		std::ofstream ofs(filename.c_str());
		ofs << ss.str();
		LOWMSG("Results written in " << filename);
	}
	if(_pTcl)
		Tcl_AppendResult(_pTcl, ss.str().c_str(), 0);
	return TCL_OK;
}

}
