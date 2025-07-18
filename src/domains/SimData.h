/*
 * SimData.h
 *
 *  Created on: Feb 15, 2011
 *      Author: ignacio.martin@imdea.org
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

#ifndef SIMDATA_H_
#define SIMDATA_H_

#include "kernel/Material.h"
#include "kernel/ParticleType.h"
#include "domains/MeshElementIterator.h"
#include "io/OutData.h"
#include "okmc/Cluster.h"
#include "SimDataAux.h"

#include <map>

namespace IO { class FileParameters; }

namespace Domains {

class SimData {

public:
	SimData(IO::FileParameters *);
	~SimData() {}

	void getACInterfaceMean(Kernel::Coordinates &, Kernel::Coordinates &, Kernel::Coordinates &) const;
	void getACInterfaceStdev(Kernel::Coordinates &, Kernel::Coordinates &, Kernel::Coordinates &) const;
	void getACInterfaceMinMax(Kernel::Coordinates &, Kernel::Coordinates &, Kernel::Coordinates &, Kernel::Coordinates &) const;
	float getACInterfaceCoverage(Kernel::Coordinates &, Kernel::Coordinates &) const;

	std::string getHistogram(Kernel::M_TYPE orig_mt, const std::string & defs) const;
	unsigned count(Kernel::M_TYPE orig_mt, const std::string &part, const std::string &defect,
			const std::string &sID, unsigned minSize, bool bDefect, const std::string &st) const;
	unsigned countPos(Kernel::M_TYPE orig_mt, const std::string &spos,
			const std::string &defe, const std::string &sID) const;
	double   radius(Kernel::M_TYPE orig_mt, const std::string &defect, const std::string &ID, double minRad) const;
	void     collectDefects(std::map<std::string, unsigned> maps[Kernel::MAX_MATERIALS]) const;
	double getDiffusivity(const std::string &mate, const std::string &name, bool bMacro) const;
        double getDose(Kernel::M_TYPE mate);
        
	IO::OutDataVectorC<double> getJumps(const std::string &name) const;
	IO::OutDataVectorC<double> getProfileMobile(const std::string &name, const std::string &st, const std::string &mate) const;

private:
	void gatherSpecificAtomFromInterfaceAndCluster(IO::OutDataVectorP<double> &odvp, Domains::MeshElementIterator const& it,
		Kernel::M_TYPE mt, Kernel::P_TYPE pt) const;
	void gatherSpecificAtomFromSpecificClusterFamily(IO::OutDataVectorP<double> &odvp, Domains::MeshElementIterator const& it,
		Kernel::M_TYPE mt, Kernel::P_TYPE pt, std::string const& def) const;
    void gatherClusterFamily(IO::OutDataVectorP<double> &odvp, Domains::MeshElementIterator const& it, Kernel::M_TYPE mt,
		std::string const& def, std::set<OKMC::Cluster const*> &was) const;
	void gatherSpecificClusterFromClusterFamily(IO::OutDataVectorP<double> &odvp, Domains::MeshElementIterator const& it, 
		Kernel::M_TYPE mt, std::string const& name, std::string const& def, std::set<OKMC::Cluster const*> &was) const;
	void gatherVanillaParticle(IO::OutDataVectorP<double> &odvp, Domains::MeshElementIterator const& it, 
		Kernel::M_TYPE mt, Kernel::P_TYPE pt, std::string const& name, std::string const& st, std::string const& def) const;
	void gatherAllInActive(IO::OutDataVectorP<double> &odvp, Domains::MeshElementIterator const& it, Kernel::P_TYPE pt) const;

public:
	IO::OutDataVectorC<double> getParticleProfile(const std::string &name, const std::string &def,
			const std::string &st, const std::string &mate) const;
	IO::OutDataVectorC<double> getAtomProfile(const std::string &name, const std::string &mate) const;
	IO::OutDataVectorC<double> getStrain(const std::string &name) const;
	IO::OutDataVectorC<double> getStress(const std::string &name) const;
        IO::OutDataVectorC<double> getBalance(const std::string &name, const std::string &mate) const;


	IO::OutDataVectorC<double> getElectrostatic(const std::string &name) const;

	std::vector<float> getMixEnthalpy(const std::string &name) const;

	std::string getCoordination(double r, const Kernel::Coordinates &c, int typ) const;
	std::string getFuzz() const; //returns an y,z surface map.

	void getLKMCInterface(std::vector<ParticleData> &c, bool bDefective) const; //LKMC interface
	void getOKMCParticles(std::vector<ParticleData> &c, const std::vector<std::string> &) const; //OKMC particles

	IO::OutDataVectorC<double> getLKMCProfile(const std::string &) const;

	IO::OutDataVectorC<double> amorphousFraction() const;
	IO::OutDataVectorC<double> getDamageProfile() const;

	double getDeltaTime() const;
	void reset();
	void setDimension(unsigned d) { _dim = d; }
	void setMinMax(const Kernel::Coordinates &m, const Kernel::Coordinates &M) { _min=m; _max=M; }
	unsigned countDisplaced(const std::string &mate, const bool bAmorph, const float threshold) const;

	static void createDefectList(Kernel::M_TYPE, const std::string &, std::vector<unsigned> &); //transform a string into list of defects
private:
	IO::OutDataVectorC<double> collapse(const IO::OutDataVectorP<double> &odv3) const;
	IO::OutDataVectorC<double> collapse(const IO::OutDataVectorC<double> &odv3) const;
	unsigned toODVC(float d[3], unsigned n[3], const Kernel::Coordinates &c) const;

	double   _lastTime;
	unsigned _dim;
	float    _odelta[3];
	Kernel::Coordinates _min, _max; //user option to filter outputs to this region.
};

}

#endif /* SIMDATA_H_ */
