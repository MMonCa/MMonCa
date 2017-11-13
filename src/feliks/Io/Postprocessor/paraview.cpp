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
/* paraviewinterface.c
 *
 * i. romero, nov 2008
 *
 * the functions that enable the communication of data with the post processor paraview
 *
 * Notes:
 
 1) It is important to determine the endian format when employing binary data. This is done in the feliks.cpp
 
 2) When writing binary data encoded in 64bits, it is necessary to put all the data in an array,
 encode it, and then dump it. It is wrong to encode each double, and write them one after another.
 Only if the number of doubles is a multiple of three the two processes give the same answer.
 */


#include "paraview.h"

#include <vector>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>

#include "Io/usercommand.h"
#include "boost/foreach.hpp"

#include "Io/base64.h"
#include "Elements/element.h"
#include "Elements/evalspot.h"
#include "General/feliksutil.h"
#include "Math/Topology/topology.h"
#include "Io/message.h"
#include "Io/logger.h"

#include "Model/model.h"
#include "Model/Parts/modelpart.h"
#include "Model/Parts/body.h"
#include "Model/Parts/poorbody.h"
#include "Model/Node/node.h"
#include "smoother.h"
#include "Math/tensor.h"


#ifdef WITHTBB
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#endif


using namespace std;
using namespace blue;


#ifdef WITHTBB
// data structure for tbb threaded version
class tbbDumpProyection{
public:
	tbbDumpProyection(std::vector<node*>&   nds, 
                      smoother              *theSmootherc, 
                      const resultdata&     r,
                      body&                 thePart) 
	: 
    my_nodes(&nds),
    theSmoother(theSmootherc), 
    grad(r),
    my_part(&thePart)
    {}
	
	void operator()( const tbb::blocked_range<size_t>& r) const
	{
        BOOST_FOREACH(node* p, *my_nodes)
		{
			p->theSmoothedResult = grad ; 
			theSmoother->smoothToNode(p->theSmoothedResult, *p, *my_part);
		}
	}
	
private:
    std::vector<node*> *const   my_nodes;
	smoother*                   theSmoother;
	resultdata                  grad;
    body*                       my_part;
};
#endif



paraview :: paraview() :
isMeshFileInitialized(false),
isResultFileInitialized(false),
areNodalCoordinatesDumped(false),
areGaussPointsInitialized(false),
step(0),
nodecounter(0),
buffer(0),
buffersize(0),
endianFormat("")
{
	name = "Data output with format for Paraview postprocessor.";
	
	// endian format
	extern bool global_bigEndianSystem;
	if (global_bigEndianSystem) endianFormat = "BigEndian";
	else						endianFormat = "LittleEndian";
	
	
	mainfile.open("feliks.pvd");
	if (!mainfile) ErrorMessage("Error opening paraview file");
	
	mainfile << "<?xml version=\"1.0\"?>" 
	<< "\n<VTKFile type=\"Collection\" version=\"0.1\" " << "byte_order=\"" << endianFormat << "\">"
	<< "\n<Collection>\n";
	
	//if ( directoryExists("felikspost.d") )
	std::ifstream check("felikspost.d");
	if (!check) system("mkdir felikspost.d");
	else		system("rm -f felikspost.d/*");
}



paraview :: paraview(const commandLine& cl) :
postprocessor(cl),
isMeshFileInitialized(false),
isResultFileInitialized(false),
areNodalCoordinatesDumped(false),
areGaussPointsInitialized(false),
step(0),
nodecounter(0),
buffer(0),
buffersize(0),
endianFormat("")
{
	name = "Data output with format for Paraview postprocessor.";
	
	// endian format
	extern bool global_bigEndianSystem;
	if (global_bigEndianSystem) endianFormat = "BigEndian";
	else						endianFormat = "LittleEndian";
	
	
	mainfile.open("feliks.pvd");
	if (!mainfile) ErrorMessage("Error opening paraview file");
	
	mainfile << "<?xml version=\"1.0\"?>" 
	<< "\n<VTKFile type=\"Collection\" version=\"0.1\" " << "byte_order=\"" << endianFormat << "\">"
	<< "\n<Collection>\n";
	
	//if ( directoryExists("felikspost.d") )
	std::ifstream check("felikspost.d");
	if (!check) system("mkdir felikspost.d");
	else		system("rm -f felikspost.d/*");
    
    //check if evalspots are to be printed for paraview
    if (evprint==true) {
        std::ifstream check("felikspost.evp");
        if (!check) system("mkdir felikspost.evp");
        else		system("rm -f felikspost.evp/*");
    }
}



paraview :: ~paraview()
{
	mainfile << "</Collection>"
	<< "\n</VTKFile>";
	mainfile.close();
	
	if (buffer != 0) delete [] buffer;
}



// this funcction returns true if dirname is the name of an existing directory or file
// 
bool paraview :: directoryExists(const std::string& dirname)
{
	bool ret;
	if ( access( dirname.c_str(), 0 ) == 0 )
    {
        struct stat status;
        stat( dirname.c_str(), &status );
		
        if ( status.st_mode & S_IFDIR )  ret = true;
        else
        {
			// the dirname is really a file
			ret = true;
        }
    }
    else ret = false;
	return ret;
}



void paraview :: dumpDofs(const model& m, std::ostream& of)
{
	if (m.getDofsPerNode() == 1)
	{
		of  << "\n                <DataArray type=\"Float64\" Name=\"Thermo_variable\" NumberOfComponents=\"" << 1 << "\"";
		
		// if the data is compressed and inline, we must add a header (not documented) with the blocksize
		// each double costs 8bytes. 
		// we must dump the whole array of temperatures as if it were a simple array of size 'nodecounter'
		
		if (compress)
		{
			of << " format=\"binary\">\n";
			
			// we compute the length of the whole encoded string, which we will dump in blocks of three doubles
			// at a time this is equal to 24 bytes, equal to 32 characters in the encoded string
			uint32_t length = base64::encoded_size( nodecounter * sizeof(double) );		
			std::string elength ( base64::encode(reinterpret_cast<uint8_t*>(&length), sizeof(uint32_t) ) );
			of.write(elength.data(), elength.size() );
		}
		else
			of << " format=\"ascii\">";
		
		int k=0;
		double buffer[3];
        BOOST_FOREACH(modelpart* bd, m.theParts )
        {
            BOOST_FOREACH(node* nnd, bd->nodes)
            {
                Pnode& nd=dynamic_cast<Pnode&>(*nnd);
				
                double Un1 = nd.getPDS().displacement(dofset::tn1);
                if (compress)
				{
					// we store the pressures in a buffer of size 3
					buffer[ k%3] = Un1;
					k++;
					if ( k%3 == 0 ) 
					{
						// when the buffer is full, we encode it and we dump it
						std::string encoded ( base64::encode(reinterpret_cast<uint8_t*>(buffer), 3*sizeof(double) ) );
						of.write(encoded.data(), encoded.size() ) ;	
					}
				}
				else
					of << "\n\t\t\t " << setw(9) << std::showpos << std::scientific << std::showpoint << Un1 ;
                
			}
		}
		
		// take care of the semi-full buffers
		if (compress && k%3 != 0)
		{
			std::string encoded ( base64::encode(reinterpret_cast<uint8_t*>(buffer), (k%3)*sizeof(double) ) );
			of.write(encoded.data(), encoded.size() ) ;				
		}
		
		
		of  << "\n                </DataArray>";
	}
	
	if ( m.getDofsPerNode() >= 3)
	{
		of  << "\n                <DataArray Name=\"Displacement\" type=\"Float64\" NumberOfComponents=\"" << 3 << "\"";
		
		// if the data is compressed and inline, we must add a header (not documented) with the blocksize
		// each double costs 8bytes. If dim=3, the displacement takes 24bytes. When reencoded in base 64, 3 bytes are converted into 4 bytes
		// and thus each node ends up taking 32bytes
		// In 2d, the displacement of a node takes 16 bytes and ends up taking 24 bytes.
		// the formula is ceil(bytes/3)*4
		if (compress)
		{
			of << " format=\"binary\">\n";
			uint32_t length( nodecounter * base64::encoded_size(3*sizeof(double)) );
			
			std::string elength ( base64::encode(reinterpret_cast<uint8_t*>(&length), sizeof(uint32_t) ) );
			of.write(elength.data(), elength.size() );
		}
		else
			of << " format=\"ascii\">";
		
        BOOST_FOREACH(modelpart* bd, m.theParts )
        {
            double* Un1=0;
            double  zero[] = {0.0,0.0,0.0};
            
            BOOST_FOREACH(node* nnd, bd->nodes)
            {
				nodetypeT ntype = nnd->getNodetype();
				if (ntype == _EMPTYnode) 
					Un1 = zero;
				else
				{
					Unode* nd = dynamic_cast<Unode*>(nnd);
					Un1       = nd->getUDS().displacement(dofset::tn1).components();
				}
				
				
				if (compress)
				{
					std::string encoded ( base64::encode(reinterpret_cast<uint8_t*>(Un1), 3*sizeof(double) ) );
					of.write(encoded.data(), encoded.size() ) ;
				}
				else
				{
					of << "\n\t\t\t";
					for (int i=0; i<3; i++) of << setw(9) << std::showpos << std::scientific << std::showpoint << Un1[i] << " ";
				}
			}
		}
		of  << "\n                </DataArray>";
	}
	
	// temperature if available
	if (m.getDofsPerNode() == 4)
	{
		of  << "\n                <DataArray type=\"Float64\" Name=\"Thermo_variable\" NumberOfComponents=\"" << 1 << "\"";
		
		// if the data is compressed and inline, we must add a header (not documented) with the blocksize
		// each double costs 8bytes. 
		// we must dump the whole array of temperatures as if it were a simple array of size 'nodecounter'
		
		if (compress)
		{
			of << " format=\"binary\">\n";
			
			// we compute the length of the whole encoded string, which we will dump in blocks of three doubles
			// at a time this is equal to 24 bytes, equal to 32 characters in the encoded string
			uint32_t length = base64::encoded_size( nodecounter * sizeof(double) );		
			std::string elength ( base64::encode(reinterpret_cast<uint8_t*>(&length), sizeof(uint32_t) ) );
			of.write(elength.data(), elength.size() );
		}
		else
			of << " format=\"ascii\">";
		
		int k=0;
		double buffer[3];
        BOOST_FOREACH(modelpart* bd, m.theParts )
        {
            BOOST_FOREACH(node* nnd, bd->nodes)
            {                
				UPnode& nd=dynamic_cast<UPnode&>(*nnd);
				
				double Un1 = nd.getPDS().displacement(dofset::tn1);
				if (compress)
				{
					// we store the pressures in a buffer of size 3
					buffer[ k%3] = Un1;
					k++;
					if ( k%3 == 0 ) 
					{
						// when the buffer is full, we encode it and we dump it
						std::string encoded ( base64::encode(reinterpret_cast<uint8_t*>(buffer), 3*sizeof(double) ) );
						of.write(encoded.data(), encoded.size() ) ;	
					}
				}
				else
				{
					of << "\n\t\t\t";
					of << " " << setw(9) << std::scientific << std::showpoint << Un1 ;
				}
			}
		}
		
		// take care of the semi-full buffers
		if (compress && k%3 != 0)
		{
			std::string encoded ( base64::encode(reinterpret_cast<uint8_t*>(buffer), (k%3)*sizeof(double) ) );
			of.write(encoded.data(), encoded.size() ) ;				
		}
		
		
		of  << "\n                </DataArray>";
	}
}



void paraview :: dumpEvpDisplacements(const model& m, std::ostream& of)
{
    
    // write header describing all nodal data. 
	// Displacements is a must
	of  << "\n            <PointData Vectors=\"Displacement\" >";
	
	if ( m.getDofsPerNode() >= 3)
	{
		of  << "\n                <DataArray Name=\"Displacement\" type=\"Float64\" NumberOfComponents=\"" << 3 << "\"";
		
		// if the data is compressed and inline, we must add a header (not documented) with the blocksize
		// each double costs 8bytes. If dim=3, the displacement takes 24bytes. When reencoded in base 64, 3 bytes are converted into 4 bytes
		// and thus each node ends up taking 32bytes
		// In 2d, the displacement of a node takes 16 bytes and ends up taking 24 bytes.
		// the formula is ceil(bytes/3)*4
		if (compress)
		{
			of << " format=\"binary\">\n";
			uint32_t length( nodecounter * base64::encoded_size(3*sizeof(double)) );
			
			std::string elength ( base64::encode(reinterpret_cast<uint8_t*>(&length), sizeof(uint32_t) ) );
			of.write(elength.data(), elength.size() );
		}
		else
			of << " format=\"ascii\">";
		
        BOOST_FOREACH(modelpart* bd, m.theParts )
        {
            double* Un1=0;
            double  zero[] = {0.0,0.0,0.0};
            
            BOOST_FOREACH(evalspot* evp, bd->evalspots)
            {
				Un1 = zero;
                ivector tmp=evp->coordinates()-evp->getReferenceCoordinates();
				Un1 = tmp.components();
                
				
				if (compress)
				{
					std::string encoded ( base64::encode(reinterpret_cast<uint8_t*>(Un1), 3*sizeof(double) ) );
					of.write(encoded.data(), encoded.size() ) ;
				}
				else
				{
					of << "\n\t\t\t";
					for (int i=0; i<3; i++) of << setw(9) << std::showpos << std::scientific << std::showpoint << Un1[i] << " ";
				}
			}
		}
		of  << "\n                </DataArray>";
	}
	
	of  << "\n            </PointData>";	
    
}

void paraview :: dumpAllElements(const model& m, std::ostream& of)
{
	of  << "            <Cells>";
	of  << "\n                <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
	
	{
        BOOST_FOREACH(modelpart* bd, m.theParts )
        {
            BOOST_FOREACH(element* ele, bd->elements)
            {
				element &el(*ele);
				
				if (el.getGeometryType() != NOPLOT)
				{
					of << "\t\t\t";
					for (int i=0; i< el.getNNodes() ; i++) 
                        of << " " << setw(6) << noshowpos << labelToPosition[ el.getNode(i).getLabel() ];
					of << "\n";
				}
			}
		}
	}
	
	
	of  << "                </DataArray>\n";
	
	int count(0);
	of  << "                <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	of << "\t\t\t";
	{
        BOOST_FOREACH(modelpart* bd, m.theParts )
        {
            BOOST_FOREACH(element* ele, bd->elements)
            {
				element &el(*ele);
				if (el.getGeometryType() != NOPLOT)
				{
					count += el.getNNodes();
					of << " " << count;
				}
			}
		}
	}
    
	of  << "\n                </DataArray>";
	
	of  << "\n                <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
	of  << "\n\t\t\t";
	{
        BOOST_FOREACH(modelpart* bd, m.theParts )
        {
            BOOST_FOREACH(element* ele, bd->elements)
            {
				element &el(*ele);
				if (el.getGeometryType() != NOPLOT) of << " " << vtkElementType(el);
			}
		}
	}
    
	of  << "\n                </DataArray>";
	of  << "\n            </Cells>";
	of  << "\n            <CellData>";
	of  << "\n                <DataArray type=\"UInt8\" Name=\"Elmt types\" format=\"ascii\">\n";
	
    {
        BOOST_FOREACH(modelpart* bd, m.theParts )
        {
            BOOST_FOREACH(element* ele, bd->elements)
            {
				element &el(*ele);				
				if (el.getGeometryType() != NOPLOT) of <<  " " << el.getEltype().label();
			}
		}
	}
	
    
	of  << "\n                </DataArray>\n";
	of  << "                <DataArray type=\"UInt8\" Name=\"Material\" format=\"ascii\">\n";
	{
        BOOST_FOREACH(modelpart* bd, m.theParts )
        {
            BOOST_FOREACH(element* ele, bd->elements)
            {
				element &el(*ele);
				if (el.getGeometryType() != NOPLOT) of <<  " " << el.getEltype().materialLabel();
			}
		}
	}
	
	of  << "\n                </DataArray>";
	of  << "\n            </CellData>";
}




void paraview :: dumpElements(const model& m, std::ostream& of)
{
    
	of  << "            <Cells>";
	of  << "\n                <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    
    int nodedisc = 0;
    
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(element* ele, bd->elements)
        {
			element &el(*ele);
			
			if (el.getGeometryType() != NOPLOT)
			{
				of << "\t\t\t" << setw(6) << noshowpos;
				for (int i=0; i< el.getNNodes() ; i++) of << " " << labelToPosition[ el.getNode(i).getLabel() ];
				of << "\n";
			}
		}
        
        if ( bd->elements.empty() ) 
        {
            // Then body treated just as points.
            int jump = 0;
            of << "\t\t\t" << setw(6) << noshowpos;
            for (int j = 0; j < bd->nodes.size(); j++) {
                of << "  "<< j +nodedisc;
                jump++;
                if ( jump >=15)
                {
                    of<<" \n";
                    of << "\t\t\t" << setw(6) << noshowpos;
                    jump=0;
                }
            }
            nodedisc += bd->nodes.size();    
        }
        
	}
    nodedisc = 0;
	of  << "                </DataArray>\n";
	
	int count(0);
	of  << "                <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	of << "\t\t\t";
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(element* ele, bd->elements)
        {
			element &el(*ele);
			if (el.getGeometryType() != NOPLOT)
			{
				count += el.getNNodes();
				of << " " << count;
			}
		}
        if ( bd->elements.empty() ) {
            // Then body treated just as points.
            int jump = 0;
            of << "\t\t\t" << setw(6) << noshowpos;
            for (int j = 0; j < bd->nodes.size(); j++) {
                of << "  "<< j+1 + nodedisc;
                jump++;
                if ( jump >=15)
                {
                    of<<" \n";
                    of << "\t\t\t" << setw(6) << noshowpos;
                    jump=0;
                }
            } 
            nodedisc += bd->nodes.size();
        }
	}
    nodedisc = 0;
	of  << "\n                </DataArray>";
	
	of  << "\n                <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
	of  << "\n\t\t\t";
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(element* ele, bd->elements)
        {
			element &el(*ele);
			if (el.getGeometryType() != NOPLOT) of << " " << vtkElementType(el);
		}
        if ( bd->elements.empty() ) {
            // Then body treated just as points.
            int jump = 0;
            of << "\t\t\t" << setw(6) << noshowpos;
            for (int j = 0; j < bd->nodes.size(); j++) {
                of << "  "<< 1 ;
                jump++;
                if ( jump >=15)
                {
                    of<<" \n";
                    of << "\t\t\t" << setw(6) << noshowpos;
                    jump=0;
                }
            }
        }
	}
    nodedisc = 0;
    
	of  << "\n                </DataArray>";
	of  << "\n            </Cells>";
	of  << "\n            <CellData>";
	of  << "\n                <DataArray type=\"UInt8\" Name=\"Elmt types\" format=\"ascii\">\n";
    
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(element* ele, bd->elements)
        {
			element &el(*ele);
			if (el.getGeometryType() != NOPLOT) of <<  " " << el.getEltype().label();
		}
        if ( bd->elements.empty() ) {
            // Then body treated just as points.
            int jump = 0;
            of << "\t\t\t" << setw(6) << noshowpos;
            for (int j = 0; j < bd->nodes.size(); j++) {
                of << "  "<< 1 + nodedisc;
                jump++;
                if ( jump >=15)
                {
                    of<<" \n";
                    of << "\t\t\t" << setw(6) << noshowpos;
                    jump=0;
                }
            } 
            nodedisc ++;
        }
	}
    nodedisc = 0;
	of  << "\n                </DataArray>\n";
	
	of  << "                <DataArray type=\"UInt8\" Name=\"Color\" format=\"ascii\">\n";
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(element* ele, bd->elements)
        {
			element &el(*ele);
			if (el.getGeometryType() != NOPLOT) of <<  " " << el.getNNodes();
		}
        if ( bd->elements.empty() ) {
            // Then body treated just as points.
            int jump = 0;
            of << "\t\t\t" << setw(6) << noshowpos;
            for (int j = 0; j < bd->nodes.size(); j++) {
                of << "  "<< 1 + nodedisc;
                jump++;
                if ( jump >=15)
                {
                    of<<" \n";
                    of << "\t\t\t" << setw(6) << noshowpos;
                    jump=0;
                }
            } 
            nodedisc++;
        }
	}
	of  << "\n                </DataArray>";
	of  << "\n            </CellData>";
}



void paraview :: dumpFakeElements(const model& m, std::ostream& of)
{
 	of  << "            <Cells>";
    
    //connectivities
	of  << "\n                <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    
    std::vector<int> evpsize_each_part;
    evpsize_each_part.push_back(0);
    
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        int lag = 0;
        for (int a = 0; a < evpsize_each_part.size(); a++) {
            lag += evpsize_each_part[a];
        }
        evpsize_each_part.push_back(bd->evalspots.size());
        
        int jump = 0;
        of << "\t\t\t" << setw(6) << noshowpos;
        for (int j = 0; j < bd->evalspots.size(); j++) 
        {
            of << "  "<< j + lag;
            jump++;
            if ( jump >=15)
            {
                of<<" \n";
                of << "\t\t\t" << setw(6) << noshowpos;
                jump=0;
            }
        }
    }
    
	of  << "                </DataArray>\n";
	
    evpsize_each_part.clear();
    evpsize_each_part.push_back(0);
    
    
    //offsets
	of  << "                <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
	of << "\t\t\t";
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        int lag = 0;
        for (int a = 0; a < evpsize_each_part.size(); a++) {
            lag += evpsize_each_part[a];
        }
        evpsize_each_part.push_back(bd->evalspots.size());
        
        int jump = 0;
        of << "\t\t\t" << setw(6) << noshowpos;
        for (int j = 0; j < bd->evalspots.size(); j++) 
        {
            of << "  "<< j+lag+1;
            jump++;
            if ( jump >=15)
            {
                of<<" \n";
                of << "\t\t\t" << setw(6) << noshowpos;
                jump=0;
            }
        } 
	}
	of  << "\n                </DataArray>";
	
    
    //types
	of  << "\n                <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
	of  << "\n\t\t\t";
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        int jump = 0;
        of << "\t\t\t" << setw(6) << noshowpos;
        for (int j = 0; j < bd->evalspots.size(); j++) {
            of << "  "<< 1;
            jump++;
            if ( jump >=15)
            {
                of<<" \n";
                of << "\t\t\t" << setw(6) << noshowpos;
                jump=0;
            }
        } 
    }
	of  << "\n                </DataArray>";
	of  << "\n            </Cells>";
	
    
    of  << "\n            <CellData>";
    of  << "\n                <DataArray type=\"UInt8\" Name=\"Elmt types\" format=\"ascii\">\n";
    
    int nodedisc = 0;
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        if ( bd->elements.empty() ) {
            // Then body treated just as points.
            int jump = 0;
            of << "\t\t\t" << setw(6) << noshowpos;
            for (int j = 0; j < bd->evalspots.size(); j++) {
                of << "  "<< 1 + nodedisc;
                jump++;
                if ( jump >=15)
                {
                    of<<" \n";
                    of << "\t\t\t" << setw(6) << noshowpos;
                    jump=0;
                }
            } 
            nodedisc ++;
        }
    }
    nodedisc = 0;
    of  << "\n                </DataArray>\n";
    
    of  << "                <DataArray type=\"UInt8\" Name=\"Color\" format=\"ascii\">\n";
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        if ( bd->elements.empty() ) {  // Then body treated just as points.
            int jump = 0;
            of << "\t\t\t" << setw(6) << noshowpos;
            for (int j = 0; j < bd->evalspots.size(); j++) {
                of << "  "<< 1 + nodedisc;
                jump++;
                if ( jump >=15)
                {
                    of<<" \n";
                    of << "\t\t\t" << setw(6) << noshowpos;
                    jump=0;
                }
            } 
            nodedisc++;
        }
    }
    of  << "\n                </DataArray>";
    of  << "\n            </CellData>";
    
}


bool paraview :: dumpMesh(const model &m) 
{
	return false;
}




bool paraview :: dump(const feliks::meshing::solidmesh &m) const
{
    return false;
}


void paraview :: dumpEvpRefCoord(const model& m, std::ostream& of)
{
	int ndm(3);
	labelToPosition.clear();
	nodecounter = 0;
	
	of  << "            <Points>";
	of  << "\n                <DataArray type=\"Float64\" NumberOfComponents=\"" << ndm << "\"";
	
	
	if (compress)
	{
		of << " format=\"binary\">\n";
		uint32_t    length(  base64::encoded_size(  m.getNEvalspots() * sizeof(double)* ndm) );
		std::string elength ( base64::encode(reinterpret_cast<uint8_t*>(&length), sizeof(uint32_t) ) );
		of.write(elength.data(), elength.size() );
	}
	else
		of << " format=\"ascii\">\n";
	
	if ( compress )
	{
		if ( buffer == 0)
		{
			buffersize = m.getNEvalspots()*ndm;
			buffer = new double[buffersize];
		}
		else if (buffersize < m.getNEvalspots()*ndm)
		{
			delete [] buffer;
			buffersize = m.getNEvalspots()*ndm;
			buffer = new double[buffersize];
		}
	}	
	
	
	int k(0);
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(evalspot* evp, bd->evalspots)
        {  
            const ivector& coor= evp->getReferenceCoordinates();
			if (compress)
			{
				for (int a=0; a<ndm; a++) buffer[k++] = coor(a);
			}
			else 
			{
				of << "\t\t\t";
				for (int i=0; i<ndm; i++) 
                    of << " " << setw(9) << std::showpos << std::scientific << std::showpoint << coor(i);
				of << "\n";
			}
            
		}
	}
	
	if (compress)
	{
		std::string encoded ( base64::encode(reinterpret_cast<uint8_t*>(buffer), sizeof(double) * m.getNEvalspots() * ndm ) );
		of.write(encoded.data(), encoded.size() ) ;
	}
	
	of  << "                </DataArray>\n";
	of  << "            </Points>\n";
}


/* Dump vibrations modes. If there are elements whose nodes have more than 4 dofs,
 it is supossed that these elements are 3D beams. 
 */
bool paraview :: dumpModes(const model &m, matrix evectors)
{
	return true;
}





bool paraview :: dumpNodalData(const model &m, std::ostream& of) 
{	
	// write header describing all nodal data. 
	// Displacements is a must
	of  << "\n            <PointData Vectors=\"Displacement\"";
	
	if ( m.getDofsPerNode() == 4) of << " Scalar=\"Thermo_variable\"";
	of << ">";
	dumpDofs(m, of);
	
	if ( results.empty()) 
	{
		of  << "\n            </PointData>";	
		return false;
	}
	
	
	// Then go possible rates
	std::vector<std::string>::iterator iter = results.begin();
	while (iter != results.end())
	{ 
		if		( *iter == "velocity" )			dumpRates(m, of, "velocity");
		else if ( *iter == "acceleration") 		dumpRates(m, of, "acceleration");
		++iter;
	}
	
	// Then go projected quantities
	iter = results.begin();
	while (iter != results.end() )
	{
		if		( resultdata::getResultType(*iter) == RESULT_TYPE_SCALAR ) 
			dumpProjection(m, of, *iter);
		else if ( resultdata::getResultType(*iter) == RESULT_TYPE_VECTOR )
			dumpProjection(m, of, *iter);
		else if ( resultdata::getResultType(*iter) == RESULT_TYPE_TENSOR )
			dumpProjection(m, of, *iter);
		else
			logger::mainlog << "\n Error in result data. Type unknown";
		++iter;
	}
	of  << "\n            </PointData>";	
	
	
	return true;
}





// while the nodes are being dumped, a map is created associating (label, printing position),
// because feliks allows negative and non-consecutive labels, but the reference in the elements
// must be to their printing position
// the key for searching is the nodal label, and the value is its printing position
void paraview :: dumpNodes(const model& m, std::ostream& of)
{
	int ndm(3);
	labelToPosition.clear();
	nodecounter = 0;
	
	of  << "            <Points>";
	of  << "\n                <DataArray type=\"Float64\" NumberOfComponents=\"" << ndm << "\"";
	
	
	if (compress)
	{
		of << " format=\"binary\">\n";
		uint32_t    length(  base64::encoded_size(  m.getNNodes() * sizeof(double)* ndm) );
		std::string elength ( base64::encode(reinterpret_cast<uint8_t*>(&length), sizeof(uint32_t) ) );
		of.write(elength.data(), elength.size() );
	}
	else
		of << " format=\"ascii\">\n";
	
	if ( compress )
	{
		if ( buffer == 0)
		{
			buffersize = m.getNNodes()*ndm;
			buffer = new double[buffersize];
		}
		else if (buffersize < m.getNNodes()*ndm)
		{
			delete [] buffer;
			buffersize = m.getNNodes()*ndm;
			buffer = new double[buffersize];
		}
	}	
	
	
	int k(0);
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(node* nnd, bd->nodes)
        {            
            const ivector &coor = nnd->getReferenceCoordinates();
			if (compress)
			{
				for (int a=0; a<ndm; a++) buffer[k++] = coor(a);
			}
			else 
			{
				of << "\t\t\t";
				for (int i=0; i<ndm; i++) 
                    of << " " << setw(9) << std::showpos << std::scientific << std::showpoint << coor(i);
				of << "\n";
			}
			labelToPosition.insert( std::make_pair( nnd->getLabel() , nodecounter++ ) );
		}
	}
	
	if (compress)
	{
		std::string encoded ( base64::encode(reinterpret_cast<uint8_t*>(buffer), sizeof(double) * m.getNNodes() * ndm ) );
		of.write(encoded.data(), encoded.size() ) ;
	}
	
	of  << "                </DataArray>\n";
	of  << "            </Points>\n";
}





class ParaviewDumpWorker
{
private:
    smoother*   _theSmoother;
    resultdata  _theResult;
    poorbody*   _thePoorbody;
    
public:
	ParaviewDumpWorker(smoother* theSmoother, resultdata& theResult, poorbody& theBody)
    :
    _theSmoother(theSmoother),
    _theResult(theResult),
    _thePoorbody(&theBody)
    {}
    
    
	void operator()(node* nd) const
	{
        nd->theSmoothedResult = _theResult;
        if (_thePoorbody != 0) {
            _theSmoother->smoothToNode(nd->theSmoothedResult, *nd, *_thePoorbody);
        }
	}
};




// writes in paraview format the gradients and other projected quantities (stresses, ...)
void paraview :: dumpProjection(const model& m, std::ostream& of, const std::string& name)
{
	resultdata grad(name);
	
	int ncomponents(resultdata::getResultNComponents(resultdata::getResultCode(name)) );
	
	of  << "\n                <DataArray Name=\"" << name << "\" type=\"Float64\"";
	of	<< " NumberOfComponents=\"" << std::noshowpos <<  ncomponents << "\"";
	
	if (compress)
	{
		of << " format=\"binary\">";
		uint32_t    length(  base64::encoded_size(  nodecounter * sizeof(double)* ncomponents) );
		std::string elength ( base64::encode(reinterpret_cast<uint8_t*>(&length), sizeof(uint32_t) ) );
		of.write(elength.data(), elength.size() );
	}
	else
		of << " format=\"ascii\">";
	
	// the binary data goes in a buffer which is encoded at the end
	int k(0);
	if ( compress )
	{
		if ( buffer == 0)
		{
			buffersize = nodecounter*ncomponents;
			buffer = new double[buffersize];
		}
		else if (buffersize < nodecounter*ncomponents)
		{
			delete [] buffer;
			buffersize = nodecounter*ncomponents;
			buffer = new double[buffersize];
		}
	}	
	
	
	
    BOOST_FOREACH(poorbody* aBody, m.thePoorbodies )
    {	
        
        feliks_for_each(aBody->nodes.begin(), aBody->nodes.end(), ParaviewDumpWorker(theSmoother, grad, *aBody));
        
        /*        
         #ifdef WITHTBB
         
         tbb::parallel_for(tbb::blocked_range<size_t>(0,aBody->nodes.size()), 
         tbbDumpProyection(aBody->nodes, theSmoother, grad, *aBody) );
         
         
         
         #else
         BOOST_FOREACH(node* nnd, aBody->nodes )
         {
         nnd->theSmoothedResult = grad ; 	
         theSmoother->smoothToNode( nnd->theSmoothedResult, *nnd, *aBody);
         }
         #endif
         */
        
        BOOST_FOREACH(node* nnd, aBody->nodes )
        {                    
            if (compress)
                for (int a=0; a<ncomponents; a++) buffer[k++] = nnd->theSmoothedResult[a];
            else 
                of << "\n\t\t\t" << setw(11) << std::showpos << std::scientific << std::right 
                << std::showpoint << nnd->theSmoothedResult;			
        }
		//}
	}
	

	if (compress)
	{
		std::string encoded(base64::encode(reinterpret_cast<uint8_t*>(buffer), sizeof(double)*nodecounter*ncomponents));
		of.write(encoded.data(), encoded.size() ) ;
	}
	
	of  << "\n                </DataArray>";
}




void paraview :: dumpRates(const model& m, std::ostream& of, const std::string& ratename)
{
    int ndm = 3;
	if ( ratename == "velocity")
		of  << "\n                <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"" 
        << ndm << "\"";
	else
		of  << "\n                <DataArray type=\"Float64\" Name=\"Acceleration\" NumberOfComponents=\"" 
        << ndm << "\"";
	if (compress)
	{
		of << " format=\"binary\">\n";
		uint32_t length( nodecounter * base64::encoded_size(ndm*sizeof(double)) );
		std::string elength ( base64::encode(reinterpret_cast<uint8_t*>(&length), sizeof(uint32_t) ) );
		of.write(elength.data(), elength.size() );
	}
	else
		of << " format=\"ascii\">";
	
	
	
	int ratenumber = ratename == "velocity" ? 1 : 2;
    BOOST_FOREACH(modelpart* bd, m.theParts )
    {
        BOOST_FOREACH(node* nnd, bd->nodes)
        {
			ivector& vrate = (ratenumber == 1) ? nnd->velocity(dofset::tn1) : nnd->acceleration(dofset::tn1);
			double *rate = vrate.components();
			
			if (compress)
			{
				std::string encoded( base64::encode(reinterpret_cast<uint8_t*>(rate), ndm*sizeof(double) ) );
				of.write(encoded.data(), encoded.size() ) ;
			}
			else
			{
				of << "\n\t\t\t";
				for (int i=0; i<ndm; i++) of << " " << setw(9) << std::scientific << std::showpoint << rate[i] ;
			}			
		}
	}
	of  << "\n                </DataArray>";
}




bool paraview :: dumpSolution2(const model &m, const double time)
{
	// write index in main file
	extern double global_tn1;
	mainfile << "   <DataSet timestep=\""<< setfill('0') << setw(5) << global_tn1  //<< step
    << "\" group=\"\" part=\"0\" file=\"felikspost.d/feliks_T"
    << setfill('0') << setw(5) << step << ".vtu\"/>\n";
	
	// write data file in vtk format
	std::stringstream filename;
	filename << "felikspost.d/feliks_T" << setfill('0') << setw(5) << step << ".vtu";
    
	std::ofstream of(filename.str().c_str(), std::ios::binary | std::ios::app);
	if (!of) std::cout << "Error opening vtu file " << filename.str() << std::endl;
	
	of  << "<?xml version=\"1.0\"?>\n";
	of  << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianFormat << "\""
    << ">" 	<< "\n    <UnstructuredGrid>";	
	
    if ( m.getNVisibleElements() != 0) {
        of  << "\n        <Piece NumberOfPoints=\"" << m.getNNodes() << "\" NumberOfCells=\"" << m.getNVisibleElements() << "\">\n";
    }
	else
    {
        of  << "\n        <Piece NumberOfPoints=\"" << m.getNNodes() << "\" NumberOfCells=\"" << m.getNNodes() << "\">\n";   
    }
	//write nodes
	dumpNodes(m, of);
	
	//write elements
	dumpElements(m, of);
	
	//write elements
	//dumpAllElements(m, of); //Including Contact Elements
	
	
	//write solution
	dumpNodalData(m, of);
	
	of  << "\n        </Piece>" 
	<< "\n    </UnstructuredGrid>" 
	<< "\n</VTKFile>\n";
	
    if (evprint == true) {
        dumpEvp2(m,time);
    }
    
    step++;
	return true;
}


bool paraview :: dumpEvp2(const model &m, const double time)
{
	
	// write data file in vtk format
	std::stringstream filename;
	filename << "felikspost.evp/feliks_T" << setfill('0') << setw(5) << step << ".vtu";
    
	std::ofstream of(filename.str().c_str(), std::ios::binary | std::ios::app);
	if (!of) std::cout << "Error opening vtu file " << filename.str() << std::endl;
	
	of  << "<?xml version=\"1.0\"?>\n";
	of  << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"" << endianFormat << "\""
    << ">" 	<< "\n    <UnstructuredGrid>";	
	
    if ( m.getNVisibleElements() != 0) {
        of  << "\n        <Piece NumberOfPoints=\"" << m.getNEvalspots() << "\" NumberOfCells=\"" << m.getNEvalspots() << "\">\n";
    }
	else
    {
        of  << "\n        <Piece NumberOfPoints=\"" << m.getNEvalspots() << "\" NumberOfCells=\"" << m.getNEvalspots() << "\">\n";   
    }
	//write evalspot original positions
	dumpEvpRefCoord(m, of);
	
	//write elements
	dumpFakeElements(m, of);
	
	//write elements
	//dumpAllElements(m, of); //Including Contact Elements
	
	
	//write solution
	dumpEvpDisplacements(m, of);
	
	of  << "\n        </Piece>" 
	<< "\n    </UnstructuredGrid>" 
	<< "\n</VTKFile>\n";
	
	return true;
}




void paraview :: print(ostream &of)
{
	postprocessor::print(of);
}




int paraview :: vtkElementType(const element& e)
{
	int       ret(0);
	geometryT geo   = e.getGeometryType();
	int       nnode = e.getNNodes();
	
	switch (geo)
	{
		case (NOPLOT) :
			ret = -1;
			break; // this elements are not plotted
			
		case (POINT):
			ret = 1;
			break;
			
		case (CURVE):
			if		(nnode ==2)     ret = 3;
			else if	(nnode ==3)     ret = 21;
			else                    ret = 4;
			break;
			
		case (SURFACE):
			if      (nnode == 3)	ret = 5;
			else if (nnode == 4)    ret = 9;
			else if (nnode == 6)    ret = 22;
			else if (nnode == 8)    ret = 23;
			break;
			
		case (VOLUME):
			if      (nnode == 4)    ret = 10;
			else if (nnode == 8)    ret = 12;
			else if (nnode ==10)    ret = 24;
			else if (nnode ==20)    ret = 25;
			break;
			
		default:
			ret = 0;
	}
	return ret;
}
