/*
 * Global.h
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

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <string>
#include <set>
#include <map>
#include <fstream>
#include <vector>
#include "kernel/Material.h"
#include "kernel/Coordinates.h"
#include "MeshElementIterator.h"
#include "DefectIterator.h"
#include "ParticleIterator.h"

namespace Kernel { class Domain; class UpdateManager; }
namespace IO { class FileParameters; class ParameterManager;}

struct Tcl_Interp;

namespace Domains {

class SimData;
class Splitter;
class MCClient;

class Global
{
public:
	Global();
	~Global();

	static Global * instance();
	static void out(int, const std::string &str);
	static void warning(const std::string &str);
	static void error(const std::string &str);
	static void setTcl(Tcl_Interp *p) { _pGlobalTcl = p; }
	static Tcl_Interp * getTcl() {return _pGlobalTcl; }
	static void setInputFileName(const std::string &s) { _inputFileName = s; }
	static void setLogFileName(const std::string &s) { _ofsLog.open(s.c_str()); }

	void deleteGlobal();

	bool & debug() { return _bDebug; }
	void setVerbosity(int verb) { _verbosity = verb; }
	const std::string & getInputFileName() const { return _inputFileName; }

	void setInteraction(Kernel::M_TYPE mt, const std::string &key) { _interactionUsed[mt].insert(key); }
	void checkInteractions(const IO::ParameterManager *pPM, const IO::FileParameters *pPar) const;
	bool isReaction(Kernel::M_TYPE, std::string one, std::string two) const;

	IO::FileParameters * getFileParameters() const { return _pPar; }
	void setFileParameters(IO::FileParameters *);
	IO::ParameterManager * PM() { return _pPM; }
	MCClient * client();
	SimData * data();

	Kernel::Domain * getDomain(const Kernel::Coordinates &);
	Kernel::Domain * getDomain(unsigned i) { return _domains[i]; }
	unsigned getDomains() const { return _domains.size(); }
   	const Domains::Splitter * getSplitter() const { return _pSplitter; }

	//iterators
	MeshElementIterator beginMEI();
	MeshElementIterator endMEI();
	DefectIterator      beginDI();
	DefectIterator      endDI();
	ParticleIterator    beginPI();
	ParticleIterator    endPI();

	void setClient(MCClient *p, const Kernel::Coordinates &m, const Kernel::Coordinates &);

	void  setTempK(float K);
	float getTempK() const { return _kelvin; }
	void anneal(double time, bool bDepth, float depth, long unsigned events);
	double getTime() const;
	void setTime(double);
	void printOutput(float percent);
	void eventReport() const;
	void reactionReport() const;
	void reactionInterfaceReport() const;
	void defectReport() const;
	void meshReport() const;
	void domainReport() const;

	bool IsPoisson() const { return _isPoisson; }
	std::string mechModel() const { return _mechModel; }
   	// debug functions for interfaces
   	void printInterfaces();
   	bool checkInterfaces();

   	long unsigned getNextEvents() const;
private:
	static Global * _pGlobal;
	static int _verbosity;
	static std::ofstream _ofsLog; //logFile
	static std::map<std::string, unsigned> _warnings;
	static Tcl_Interp *_pGlobalTcl;
	static std::string _inputFileName;
	Kernel::UpdateManager * _pTimeUpdate;

	IO::FileParameters * _pPar;
	IO::ParameterManager *_pPM;
	MCClient *_pClient;
	SimData *_pSimData;
	Domains::Splitter *_pSplitter;
	std::vector<Kernel::Domain *> _domains;

	bool _isPoisson;
	std::string _mechModel;

	bool _bDebug;

	//time steps
	double _kelvin;
    double _ptimeUnit, _ptimeExp;
    double _time;
    long unsigned _oldEv;
    double _oldPrintTime;
  	time_t _oldTime, _initCPUTime;
  	double _annealCPUTime;

	std::set<std::string> _interactionUsed[Kernel::MAX_MATERIALS];
};

inline Global * global() { return Global::instance(); }

}

#endif /* GLOBAL_H_ */
