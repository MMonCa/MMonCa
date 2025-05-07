/*
 * Global.cpp
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

#include "Global.h"
#include "SimData.h"
#include "MCClient.h"
#include "SimpleSplitter.h"
#include "io/FileParameters.h"
#include "io/ParameterManager.h"
#include "io/Command.h"
#include "kernel/Material.h"
#include "kernel/Domain.h"
#include "kernel/SubDomain.h"
#include "kernel/Mesh.h"
#include "kernel/RateManager.h"
#include "kernel/UpdateManager.h"
#include "okmc/ClusterParam.h"
#include "okmc/InterfaceParam.h"
#include "okmc/Interface.h"
#include "okmc/Defect.h"

#ifdef _OPENMP
#	include <omp.h>
#endif
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <iomanip>

using std::string;
using std::map;
using std::stringstream;
using std::vector;

namespace Domains {

Global::Global() : _pTimeUpdate(0), _bDebug(false), _annealCPUTime(0)
{
	const char *p = std::getenv("MCPATH");
	string path = (p? string(p) : string("./config"));

	_pPar      = new IO::FileParameters(path);
	_pPM       = 0;
	_bDebug    = _pPar->getBool("MC/General/debug");
	_pClient   = 0;
	_pSimData  = 0;
	_pSplitter = 0;
	_isPoisson = false;

	_time = 0;
	std::time(&_oldTime);
	std::time(&_initCPUTime);
	_oldEv = 0;
	_ptimeUnit = 1;
	_ptimeExp = -15;
	_oldPrintTime = 0;
	_kelvin = 23 + 273.15;
}

Global::~Global()
{
	const unsigned TOTAL_SENTENCES = 15;
	int number = _initCPUTime % TOTAL_SENTENCES;
	delete _pSimData;
	while(_domains.size())
	{
		delete _domains.back();
		_domains.pop_back();
	}
	delete _pClient;
	delete _pPM;
	delete _pPar;
	delete _pSplitter;
	delete _pTimeUpdate;

	time_t newCPUTime;
	std::time(&newCPUTime);
	LOWMSG("Time spent: " << std::difftime(newCPUTime, _initCPUTime) << "s. Annealing: " << _annealCPUTime << "s.");

	for(map<string, unsigned>::iterator it=_warnings.begin(); it != _warnings.end(); ++it)
		LOWMSG(std::right << std::setw(10) << it->second << " times: " << it->first);

	static const string messages[TOTAL_SENTENCES] = {
			"Camarón que se duerme...\n                     ... se lo lleva la corriente",
			"¡Esto sí que es ...\n                ... una chirigota!",
			"¡Adiós amigos!",
			"Milagrosamente...\n              ... he ganado",
			"¡Hasta la vista!",
			"Hay 10 tipos de personas...\n                        ... las que saben binario y las que no",
			"¡Curioso... \n        ...cuando menos!",
			"No por mucho simular...\n                    ... amanece más temprano",
			"El mayor obstáculo para encontrar la felicidad...\n                                 ... es buscarla.",
			"Lo más importante en la vida es...\n                               ... ¡Vaya! ¡Se me olvidó!",
			"La pregunta no es si hay vida después de la muerte...\n                            ... ¡sino si la hay antes!",
			"Die Welt ist ein Irrenhaus...\n                          ... aber hier ist die Zentrale!",
			"Hay todo tipo de bacterias ahí fuera...\n                                    ... esperando a fermentarme.",
			"En la vida unas veces se gana...\n                             ... y siempre se pierde.",
			"Son, observe the time...\n                     ... and fly from evil."
	};
	LOWMSG("There are " << TOTAL_SENTENCES << " total sentences. Collect them all!\n" << messages[number]);
	_ofsLog.close();
}

Global * Global::instance()
{
	if (!_pGlobal)
		_pGlobal = new Global;
	return _pGlobal;
}

void Global::deleteGlobal()
{
	delete _pGlobal;
	_pGlobal = 0;
}

void Global::setFileParameters(IO::FileParameters *pFP)
{
	delete _pPar;
	_pPar = pFP;
}

MCClient * Global::client()
{
	if(_pClient)
		return _pClient;
	ERRORMSG("Current command needs a previous 'init' command to work properly (client not build)");
	return 0;
}

SimData * Global::data()
{
	if(_pSimData)
		return _pSimData;
	ERRORMSG("Current command needs a previous 'init' command to work properly (data not build)");
	return 0;
}

void Global::setClient(MCClient *p, const Kernel::Coordinates &m, const Kernel::Coordinates &M,
           std::vector<float> const * const aLinesX, std::vector<float> const * const aLinesY, std::vector<float> const * const aLinesZ)
{
	delete _pPM;
	_pPM = new IO::ParameterManager(_pGlobalTcl, _pPar);
	_bDebug = _pPar->getBool("MC/General/debug");
	if (_pClient || _pSimData || _domains.size())
	{
		LOWMSG("Deleting old simulation");
		delete _pClient;
		delete _pSimData;
		while(_domains.size())
		{
			delete _domains.back();
			_domains.pop_back();
		}
		delete _pSplitter;
	}
	_time = 0;
	_isPoisson = _pPar->getBool("MC/Electrostatic/load.poisson");
	_mechModel = _pPar->getString("Mechanics/General/model");
	std::time(&_oldTime);
	_oldEv = 0;
	_ptimeUnit = 1;
	_ptimeExp = -15;
	_kelvin = 23 + 273.15;
	_oldPrintTime = 0;

	_pClient = p;
	// create domains using the "splitter"
        _pSplitter = new SimpleSplitter(m, M, p->getMaterial());

        if(aLinesX != nullptr && aLinesY != nullptr && aLinesZ != nullptr) {
            if(_pSplitter->getDomains() == 1u && _pSplitter->getSubDomains(0) == 1u) {
                _domains.push_back(new Kernel::Domain(0, _pGlobalTcl, p, m, M, aLinesX, aLinesY, aLinesZ));
            }
            else {
                WARNINGMSG("Lines setting in init (inhomogeneous mesh) is only allowed when domains == 1 && subdomains == 1");
            }
        }
        else {
            for(unsigned i=0; i<_pSplitter->getDomains(); ++i)
            {
                Kernel::Coordinates area, Area;
                _pSplitter->splitDomain(i, area, Area);
                _domains.push_back(new Kernel::Domain(i, _pGlobalTcl, p, area, Area));
            }
        }
	for(std::vector<Kernel::Domain *>::iterator it=_domains.begin(); it!=_domains.end(); ++it)
		(*it)->init(p->isFromStream());
	_pSimData = new SimData(_pPar);

	int nSubDom = _pSplitter->getSubDomains(0);
#ifdef _RELEASE_PARALLEL
	bool bParallel = true;
#else
	bool bParallel = false;
#endif
	if(nSubDom == 1 && bParallel)
		WARNINGMSG("This binary is compiled for parallelization, but you use only 1 thread. We suggest using the serial binary");
	if(nSubDom > 1 && !bParallel)
		WARNINGMSG("This binary is serial, and you are requesting parallelization. We suggest using the parallel binary.");
	delete _pTimeUpdate;
	_pTimeUpdate = new Kernel::UpdateManager("MC/General/snapshot");
	_pTimeUpdate->print(false);
}

Kernel::Domain * Global::getDomain(const Kernel::Coordinates &c)
{
	if(_domains.size())
		return _domains[_pSplitter->getDomain(c)];
	else
		ERRORMSG("Current command needs a previous 'init' or 'restart' command to work properly (domain not build)");
	return 0;
}

void Global::meshReport() const
{
	_domains.back()->_pMesh->print();
}

void Global::domainReport() const
{
	_domains.back()->_pMesh->printDomains();
}

void Global::defectReport() const
{
	map<string, unsigned> maps[Kernel::MAX_MATERIALS];
	Domains::global()->data()->collectDefects(maps);
	LOWMSG("----------------------------- Defect logfile --------------");
	for(Kernel::M_TYPE mt=0; mt < Domains::global()->PM()->getNMaterials(); ++mt)
	{
		if(maps[mt].empty())
			continue;
		string matName = Domains::global()->PM()->getMaterialName(mt);
		LOWMSG("---------------------------------------" << std::right << std::setw(10) << matName << " --------");
		for(map<string, unsigned>::iterator it=maps[mt].begin(); it!=maps[mt].end(); ++it)
			LOWMSG(std::left << std::setw(30) << it->first << std::right << std::setw(8) << it->second);
	}
}

void Global::eventReport() const
{
	Kernel::EventLog evLog;
	for(vector<Kernel::Domain *>::const_iterator it=_domains.begin(); it!=_domains.end(); ++it)
		for(unsigned sub=0; sub < (*it)->_pRM->getSubDomains(); ++sub)
			evLog += (*it)->_pRM->getSubDomain(sub)->_evLog;
	evLog.print();
}

void Global::reactionReport() const
{
	Kernel::ReactionLog reLog;
	for(vector<Kernel::Domain *>::const_iterator it=_domains.begin(); it!=_domains.end(); ++it)
		for(unsigned sub=0; sub < (*it)->_pRM->getSubDomains(); ++sub)
		reLog += (*it)->_pRM->getSubDomain(sub)->_reLog;
	reLog.print();
}

void Global::reactionInterfaceReport() const
{
	Kernel::ReactionLog reLog;
	for(vector<Kernel::Domain *>::const_iterator it=_domains.begin(); it!=_domains.end(); ++it)
		for(unsigned sub=0; sub < (*it)->_pRM->getSubDomains(); ++sub)
		reLog += (*it)->_pRM->getSubDomain(sub)->_reLog;
	reLog.printInterface();
}

void Global::setTempK(float K)
{
	_kelvin = K;
	for(vector<Kernel::Domain *>::iterator it=_domains.begin(); it!=_domains.end(); ++it)
		(*it)->_pRM->setTempK(_kelvin);
}

void Global::anneal(double time, bool bDepth, float depth, long unsigned event)
{
	MEDMSG("Annealing " << time);

	time_t initCPUTime;
	std::time(&initCPUTime);

	double printTime = _pTimeUpdate->getNextTime();
	double initTime  = _time;
	double endTime   = time + _time;

	if(bDepth || event > 0)
		endTime = std::numeric_limits<double>::max();
	if(event != 0 && _domains.size() > 1)
		ERRORMSG("Cannot specify number of events as stopping conditions in parallel");

	double origDepth = 0;
	for(vector<Kernel::Domain *>::iterator it = _domains.begin(); it != _domains.end(); ++it)
		origDepth += (*it)->_pRM->getDepthLA();
	origDepth /= _domains.size();
	if(origDepth > 1e35 && bDepth) //not inited...
	{
		Kernel::Coordinates m, M, mean;
		Domains::global()->client()->getCell(m, M);
		Domains::global()->data()->getACInterfaceMean(m, M, mean);
		origDepth = mean._x;
	}
	bool bStop = false;
	long unsigned nEvents = 0;
	double newDepth = 0;
	while(_time < endTime && bStop == false)
	{
		printTime = _pTimeUpdate->getNextTime();
		nEvents = 0;
		double newTime = 0;
		newDepth = 0;
		bool bFinished = true;
		double finalTime = (printTime <= endTime) ? printTime : endTime;
        #pragma omp parallel shared(bFinished) num_threads(_domains.size()) reduction(+:nEvents,newTime,newDepth)
		{
			#pragma omp for schedule(dynamic,1)
			for(unsigned nD = 0; nD < _domains.size(); ++nD)
			{
				_domains[nD]->_pRM->setTempK(_kelvin); //to update rates
				//pRM-> anneal always returns the requested time
				_domains[nD]->_pRM->anneal(finalTime - _time, bDepth, depth, event);
				if(_domains[nD]->_pRM->getSubDomain(0)->getMaxRate() != 0)
					bFinished = false;
				nEvents  += _domains[nD]->_pRM->getEvents();
				newTime  += _domains[nD]->_pRM->getTime();
				newDepth += _domains[nD]->_pRM->getDepthLA();
			}
		}
		newTime /= _domains.size();
		newDepth /= _domains.size();
		if(newDepth > 1e35)
			newDepth = origDepth;
		_time = newTime;
		if(bFinished && event == 0 && bDepth == false)
			_time = endTime;
		(*_pTimeUpdate)(_time, nEvents, newDepth);
		if(event != 0)
		{
			if(Tcl_EvalEx(_pGlobalTcl,  "snapshot", -1, 0) != TCL_OK)
				WARNINGMSG("Snapshot not defined or error.");
			printOutput(100.*float(nEvents)/event);
		}
		else if(bDepth)
		{
			if(Tcl_EvalEx(_pGlobalTcl,  "snapshot", -1, 0) != TCL_OK)
				WARNINGMSG("Snapshot not defined or error.");
			printOutput(100.*(origDepth - newDepth)/(origDepth - depth));
		}
		else if((_time - initTime)/(endTime - initTime) < 1.)
		{
			if(Tcl_EvalEx(_pGlobalTcl,  "snapshot", -1, 0) != TCL_OK)
				WARNINGMSG("Snapshot not defined or error.");
			printOutput(100.*(_time - initTime)/(endTime - initTime));
		}
		bStop = bFinished || (event > 0 && nEvents >= event) || (bDepth && newDepth < depth);
	}
	time_t endCPUTime;
	std::time(&endCPUTime);
	_annealCPUTime += std::difftime(endCPUTime, initCPUTime);
}

double Global::getTime() const
{
	return _domains.back()->_pRM->getTime();
}

void Global::setTime(double t)
{
	for(std::vector<Kernel::Domain *>::iterator it = _domains.begin(); it != _domains.end(); ++it)
		(*it)->_pRM->setTime(t);
	_time = t;
}

namespace {
//small algorithm from stackoverflow to get the total memory.
int parse(char* line)
{
	int i = strlen(line);
	while (*line < '0' || *line > '9') line++;
	line[i-3] = '\0';
	i = atoi(line);
	return i;
}

int getMemoryValue()  //in Mb
{
	FILE* file = fopen("/proc/self/status", "r");
	if(!file)
		return 0;
	int result = -1;
	char line[128];

	while (fgets(line, 128, file) != 0)
	{
		if (strncmp(line, "VmSize:", 7) == 0)
		{
			result = parse(line);
			break;
		}
	}
	fclose(file);
	return result/1024;
}
}

void Global::printOutput(float percent)
{
	time_t newTime;
	std::time(&newTime);
	double dtime = std::difftime(newTime, _oldTime);
	long unsigned nEvents = 0;
	for (vector<Kernel::Domain *>::iterator it=_domains.begin(); it!=_domains.end(); ++it)
		nEvents += (*it)->_pRM->getEvents();
	if(nEvents != _oldEv || percent == 100.)
	{
		LOWMSG2(std::right << std::fixed   << std::setw(3) << std::setprecision(0) << _kelvin-273.15 << "C ");
		LOWMSG2(std::right << std::setw(8) << std::setprecision(3) << _time << "s ");
		LOWMSG2(std::right << std::fixed   << std::setw(11) << nEvents);
		LOWMSG2(" " << std::right << std::fixed << std::setprecision(2) << std::setw(6) << percent << "% ");
		if(_time - _oldPrintTime)
			LOWMSG2(std::left << std::scientific << std::setprecision(1) << std::setw(3) << " " << (nEvents - _oldEv)/(_time - _oldPrintTime) << " s^-1 ");
		double depthLA = 0;
		for (vector<Kernel::Domain *>::iterator it=_domains.begin(); it!=_domains.end(); ++it)
			depthLA += (*it)->_pRM->getDepthLA();
		depthLA /= _domains.size();
		if(depthLA < std::numeric_limits<float>::max()*0.9/_domains.size())   // to control if setDepthLA has been ever called
			LOWMSG2(std::right << std::fixed << std::setprecision(1) << depthLA << " nm. ");
		if(dtime != 0)
			LOWMSG2(std::right << std::fixed << std::setw(7) << unsigned((nEvents-_oldEv)/dtime) << " ev/s ");
		LOWMSG(getMemoryValue() << " Mb");
		_oldEv = nEvents;
	}
	_oldTime = newTime;
	_oldPrintTime = _time;
}


void Global::out(int verbosity, const std::string &str)
{
	if(verbosity <= _verbosity)
	{
		std::cout << str;
		std::cout.flush();
		_ofsLog << str;
	}
}

void Global::warning(const std::string &str)
{
	map<string, unsigned>::iterator it = _warnings.find(str);
	if(it != _warnings.end())
		it->second++;
	else
	{
		stringstream ss;
		ss << "---------------------------- Warning -----------------------------------\n";
		ss << str << "\n";
		std::cout << ss.str();
		_ofsLog << ss.str();
		_warnings[str] = 1;
	}
}
void Global::error(const std::string &str)
{
	stringstream ss;
    ss << "FATAL error: '" << str << "'\n";
    std::cerr << ss.str();
    _ofsLog << ss.str();
    _ofsLog.close();
    std::exit(1);
}

void Global::checkInteractions(const IO::ParameterManager *pPM, const IO::FileParameters *pPar) const
{
	for(Kernel::M_TYPE mt0=0; mt0 < pPM->getNMaterials(); ++mt0)
	{
		string mat0 = pPM->getMaterialName(mt0);
		string base = mat0 + "/Models/interactions";
		IO::array<string, string> theMap = pPar->getArray(base);
		for(IO::array<string, string>::iterator it=theMap.begin(); it!=theMap.end(); ++it)
		{
			if(it->second == "false")
				continue;
			if(_interactionUsed[mt0].find(it->first) == _interactionUsed[mt0].end())
				WARNINGMSG(base << " Interaction not used: " << it->first);
		}
	}
}

bool Global::isReaction(Kernel::M_TYPE mt, string one, string two) const
{
	Kernel::Domain *pDom = _domains[0];

	vector<string> tokens_l;  //for the two incoming species
	IO::ParameterManager::getTokens(one, ':', tokens_l);
	if(tokens_l.size() == 0 || tokens_l.size() > 2)
		ERRORMSG("Unrecognized syntax in first reactant for " << one);
	unsigned defect_type_l = pDom->_pClPar->getDefectNumber(mt, tokens_l[0]);
	if(defect_type_l == pDom->_pClPar->defectSize(mt))
	{
		WARNINGMSG("Defect " << tokens_l[0] << " not recognized");
		return false;
	}

	const OKMC::EDType *edt = pDom->_pClPar->getParams(mt, defect_type_l);
	if(tokens_l.size() == 1)  // notation is ICluster+VCluster
	{
		unsigned defect_type_r = pDom->_pClPar->getDefectNumber(mt, two);
		if(defect_type_r == pDom->_pClPar->defectSize(mt)) //second defect not MC
		{
			Kernel::Event::E_TYPE ev = Domains::global()->PM()->getDefectType(two, mt);
			if(ev == Kernel::Event::INTERFACE)
			{
				Kernel::M_TYPE mt2 = global()->PM()->getMaterialNumber(two);
				vector<bool> &ins = pDom->_pIFPar->_params[std::min(mt, mt2)][std::max(mt, mt2)]->_interactWithMC;
				return ins[defect_type_l];
			}
			else
			{
				WARNINGMSG("Defect " << tokens_l[0] << " not recognized");
				return false;
			}
		}
		//second defect is MC
		return edt->_interactMC.find(defect_type_r) != edt->_interactMC.end();
	}
	else //notation is ICluster:I53+I
	{
		Kernel::ID ID = global()->PM()->getID(mt, tokens_l[1]);
		unsigned hash = edt->_hash.cluster2hash(ID);
		if(hash == edt->_hash.invalidHash())
			return false;
		Kernel::P_TYPE pt = global()->PM()->getParticleNumber(mt, two);
		if(pt == Kernel::UNDEFINED_TYPE)
			return false;
		return edt->_hash._map[hash]._interactMP[pt];
	}
}

//iterators
MeshElementIterator Global::beginMEI()
{
	MeshElementIterator mei;
	mei._nElement = 0;
	mei._nDomain = 0;
	mei._pME = getDomain(mei._nElement)->_pMesh->getElement(mei._nElement);
	return mei;
}

MeshElementIterator Global::endMEI()
{
	return MeshElementIterator();
}

ParticleIterator Global::beginPI()
{
	ParticleIterator pi;
	pi._mei = beginMEI();
	pi._pPart = pi._mei->getFirstPart();
	if(pi._pPart)
		return pi;
	++pi;
	return pi;
}

ParticleIterator Global::endPI()
{
	return ParticleIterator();
}

DefectIterator Global::beginDI()
{
	DefectIterator di(getDomain(0)->_pRM->getSubDomain(0)->begin(0));
	if(*di._pSi != getDomain(0)->_pRM->getSubDomain(0)->end(0) && dynamic_cast<OKMC::Defect *>(**di._pSi) != 0)
		return di;
	++di;
	return di;
}

DefectIterator Global::endDI()
{
	DefectIterator di(getDomain(0)->_pRM->getSubDomain(0)->begin(0));
	di._pSi.reset();
	return di;
}

void Global::printInterfaces() {
	for(Domains::DefectIterator it = Domains::global()->beginDI(); it!= Domains::global()->endDI(); ++it)
	{
		if ((*it)->getEType() == Kernel::Event::INTERFACE) {
			OKMC::Defect *d    = const_cast<OKMC::Defect *>(*it);
			OKMC::Interface *p = static_cast<OKMC::Interface *>(d);
			Kernel::M_TYPE m1 = p->getElement(0)->getMaterial();
			Kernel::M_TYPE m2 = p->getElement(1)->getMaterial();
			LOWMSG("Interface :" << Domains::global()->PM()->getMaterialName(m1) << "/" << Domains::global()->PM()->getMaterialName(m2));
		}
	}
}

bool Global::checkInterfaces() {
	bool ret = true;
	for(Domains::DefectIterator it = Domains::global()->beginDI(); it!= Domains::global()->endDI(); ++it)
	{
		if ((*it)->getEType() == Kernel::Event::INTERFACE) {
			OKMC::Defect *d    = const_cast<OKMC::Defect *>(*it);
			OKMC::Interface *p = static_cast<OKMC::Interface *>(d);
			Kernel::M_TYPE m1 = p->getElement(0)->getMaterial();
			Kernel::M_TYPE m2 = p->getElement(1)->getMaterial();
			if (m1 == m2) {
				unsigned i1, i2;
				std::string s1, s2;

				s1 = Domains::global()->PM()->getMaterialName(m1);
				s2 = Domains::global()->PM()->getMaterialName(m2);

				i1 = p->getElement(0)->getIndex();
				i2 = p->getElement(1)->getIndex();

				LOWMSG("Interface: " << s1 << "(index=" << i1 << ")/" << s2 << "(index=" << i2 << ")");
				ret = false;
			}
		}
	}
	return ret;
}

long unsigned Global::getNextEvents() const
{
	return _pTimeUpdate->getNextEvents();
}

int Global::_verbosity = 0;
Global * Global::_pGlobal = 0;
map<string, unsigned> Global::_warnings;
std::ofstream Global::_ofsLog;
string Global::_inputFileName;
Tcl_Interp * Global::_pGlobalTcl;
}
