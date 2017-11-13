#ifndef __MeshParser_PARSER__H
#define __MeshParser_PARSER__H

#include <vector>
#include <fstream>
#include "kernel/Coordinates.h"

namespace IO {

class MeshParser
{
public:
	MeshParser(const std::string &filename, float scale, float DIMX=1, float DIMY=1, float DIMZ=1);
	virtual ~MeshParser();
	std::string getMaterial(float x, float y, float z) const;
	void getCellSize(Kernel::Coordinates &m, Kernel::Coordinates &M) const;

private:
	static const unsigned HULL=20;				// Size of each side of every HULL (superposed Cartesian grid)
	const float _scale;

	struct float3 { float x; float y; float z; };
	struct int2 { int first; int second; };
	struct int3 { int first; int second; int third; };
	struct faces	{
		unsigned int code;		// MeshParser surface type
		int edge[4];			// edge index (up to 4 faces)
		int vertex[4];			// vertex index (up tp 4 vertices)
	};
	struct elements {
		unsigned int code;		// MeshParser element type
		int face[6];			// face index (up to 6 faces)
		int vertex[8];			// vertex index (up to 8 vertices)
		std::string material;	// Material
	};

	int3 FindOctree(float x, float y, float z) const;
	int  FindElement(float x, float y, float z, int3 OCT) const;
	int2 getDiagonal(int EID) const;
	float3 GetBarycentricCoords(float x, float y, float z, int EID) const;
	void Read();

	std::ifstream gridfile;
	std::string line;
	std::string end;
	std::string Material;

	unsigned nb;
	
	unsigned nb_tetrahedra;
	unsigned nb_cubes;
	unsigned nb_triangles;
	unsigned nb_regions;

	int edge_ide;
	int face_id;

	unsigned idx, idy, idz;
	float xmin,xmax,ymin,ymax,zmin,zmax; 		// Octree calculations
	float Xmin,Xmax,Ymin,Ymax,Zmin,Zmax;		// Absolute Domain size calcs
	float dimX;									// Absolute Domain dimensions
	float dimY;
	float dimZ;

	int vertex_id, edge_id, elem_id;

	float3 		Hull[2];
	int3 		Cube[2];

	float x,y,z;
	float3 rot[3];
	float Xshift, Yshift, Zshift;
	int face_code;
	int elem_code;
	unsigned nb_vertices;
	unsigned nb_edges;
	unsigned nb_faces;
	unsigned nb_elements;
	float3	* Vertices;
	int2 * Edges;
	faces *	Faces;
	elements * Elements;
	std::vector< std::vector < std::vector< std::vector <int> > > > Octree;
};

}

#endif

