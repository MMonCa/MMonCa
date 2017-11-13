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
/* paraview.h
 *
 * collect all the functions that dump data in paraview input file
 * for further postprocessing
 *
 * i. romero,  november 2008
 */

#ifndef _paraviewpost_h
#define _paraviewpost_h

#include "Math/feliksmath.h"
#include "Io/Postprocessor/postprocessor.h"

#include <ostream>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

class element;
class model;
class commandLine;
class solidmesh;

namespace feliks
{
    namespace meshing
    {
        class solidmesh;
    }
}


using namespace std;

class paraview : public postprocessor
{
public:
    paraview();
	paraview(const commandLine& cl);
	~paraview();
	
	bool			dumpMesh(const model &m);
    bool            dump(const feliks::meshing::solidmesh &m) const;
    bool			dumpModes(const model &m, matrix evectors);
    bool			dumpSolution2(const model &m, const double time);
	virtual void	print(std::ostream &of=std::cout);

	
private:
	bool	isMeshFileInitialized;
	bool	isResultFileInitialized;
	bool	areNodalCoordinatesDumped;
	bool	areGaussPointsInitialized;

	std::ofstream mainfile;
	std::map<int,int> labelToPosition;
	int		step;
	int		nodecounter;
	double* buffer;
	int		buffersize;
	std::string endianFormat;
	
	bool	directoryExists(const std::string& dirname);
    void    dumpCellResult(const model& m, std::ostream& of, const std::string& resultname);
	void	dumpDofs(const model& m, std::ostream& of);
	void	dumpElements(const model& m, std::ostream& of);
	void	dumpAllElements(const model& m, std::ostream& of);
	bool	dumpNodalData(const model &m, std::ostream& of);
	void	dumpNodes(const model& m, std::ostream& of);
	void	dumpProjection(const model& m, std::ostream& of, const std::string& name);
	void	dumpRates(const model& m, std::ostream& of, const std::string& ratename);
	int		vtkElementType(const element& e);

    //methods for optionally plotting the evalspots as sole points.
    void	dumpEvpDisplacements(const model& m, std::ostream& of);
    void    dumpEvpRefCoord(const model& m, std::ostream& of);
    bool    dumpEvp2(const model &m, const double time);
    void	dumpFakeElements(const model& m, std::ostream& of); // empty elements, to plot sole points

    
    
};


#endif


