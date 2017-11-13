#include "MeshParser.h"
#include "Diagnostic.h"
#include <iostream>
#include <cstdlib>

using namespace IO;

MeshParser::MeshParser(const std::string &filename, float scale, float DIMX, float DIMY, float DIMZ) : _scale(scale)
{
	gridfile.open(filename.c_str(), std::ifstream::in );
	if(!gridfile)
		ERRORMSG("Grid file not found: " << filename);
	dimX = DIMX;
	dimY = DIMY;
	dimZ = DIMZ;

	Read(); //some action!
}

/*MeshParser::MeshParser()
{
}*/

MeshParser::~MeshParser()
{
}

void MeshParser::getCellSize(Kernel::Coordinates &m, Kernel::Coordinates &M) const
{
	m._x = Xmin; M._x = Xmax;
	m._y = Ymin; M._y = Ymax;
	m._z = Zmin; M._z = Zmax;
}


void MeshParser::Read()
{
	LOWMSG("Processing mesh file ...");
	//std::cout << "Processing MeshParser file." << std::endl;

	nb = 0; 				// Size number initialized to 0
	gridfile >> line;

	while (line != "Vertices"){

		if (line=="nb_vertices") {
			gridfile >> line; // read '='
			gridfile >> nb_vertices;}
		if (line=="nb_edges") {
			gridfile >> line; // read '='
			gridfile >> nb_edges;}
		if (line=="nb_faces") {
			gridfile >> line; // read '='
			gridfile >> nb_faces;}
		if (line=="nb_elements") {
			gridfile >> line; // read '='
			gridfile >> nb_elements;}
		if (line=="translate") {
			gridfile >> line; gridfile >> line; // read '= ['
			gridfile >> Xshift; gridfile >> Yshift; gridfile >> Zshift;}
		if (line=="transform") {
			gridfile >> line; gridfile >> line; // read '= ['
			gridfile >> rot[0].x; gridfile >> rot[0].y; gridfile >> rot[0].z;
			gridfile >> rot[1].x; gridfile >> rot[1].y; gridfile >> rot[1].z;
			gridfile >> rot[2].x; gridfile >> rot[2].y; gridfile >> rot[2].z;
		}
		gridfile >> line;
	};

	Vertices = new float3 [nb_vertices];
	Edges = new int2 [nb_edges];
	Faces = new faces [nb_faces];
	Elements = new elements [nb_elements];


	// -------------------------------------------------------------------- 
	//			STORING VERTICES
	// --------------------------------------------------------------------

	gridfile >> line; 		// After reading 'Vertices' it should be: '(nb_vertices)'
	line.erase(line.begin());
	line.erase(line.end()-1);	
	nb = atoi(line.c_str());	// it should be now 'nb_vertices'

	if (nb != nb_vertices) { 
		ERRORMSG("Number of vertices mismatch !! : " << nb << " != " << nb_vertices);
		/*std::cout << "ERROR - Nb of vertices mismatch!!: " << nb
			  << " != " << nb_vertices << std::endl;*/
	};

	gridfile >> line; 	// it should be '{'


	// Check absolute maximum and minimum coordinates of the domain
	
	Xmin = Ymin = Zmin = 1E15;
	Xmax = Ymax = Zmax = -1E15;

	for (unsigned i=0; i<nb_vertices; i++){
		gridfile >> line; x = atof(line.c_str()) * _scale;
		gridfile >> line; y = atof(line.c_str()) * _scale;
		gridfile >> line; z = atof(line.c_str()) * _scale;

		Vertices[i].x = (rot[0].x * x + rot[0].y * y + rot[0].z * z) + Xshift * _scale;
		Vertices[i].y = (rot[1].x * x + rot[1].y * y + rot[1].z * z) + Yshift * _scale;
		Vertices[i].z = (rot[2].x * x + rot[2].y * y + rot[2].z * z) + Zshift * _scale;
	
		if (Vertices[i].x < Xmin) Xmin = Vertices[i].x;
		if (Vertices[i].x > Xmax) Xmax = Vertices[i].x;
		if (Vertices[i].y < Ymin) Ymin = Vertices[i].y;
		if (Vertices[i].y > Ymax) Ymax = Vertices[i].y;
		if (Vertices[i].z < Zmin) Zmin = Vertices[i].z;
		if (Vertices[i].z > Zmax) Zmax = Vertices[i].z;

	};
	LOWMSG(nb_vertices << "\t vertices stored.");
	//std::cout << nb_vertices << "\t vertices stored." << std::endl;

	dimX = Xmax - Xmin;	dimY = Ymax - Ymin;	dimZ = Zmax - Zmin;
	if (Xmin < 0) Xshift = -Xmin;
	else if (Xmin > 0) Xshift = -Xmin;
	else Xshift = 0;

	if (Ymin < 0) Yshift = -Ymin;
	else if(Ymin > 0) Yshift = -Ymin;
	else Yshift = 0;

	if (Zmin < 0) Zshift = -Zmin;
	else if (Zmin > 0) Zshift = -Zmin;
	else Zshift = 0;
	
	LOWMSG("\t Diagonal vertices: MIN(" << Xmin << "," << Ymin << "," << Zmin
				       << ") <--> MAX(" << Xmax << "," << Ymax << "," << Zmax
				       << ")");
	/*std::cout << "\t Diagonal vertices: MIN(" << Xmin << "," << Ymin << "," << Zmin
			      << ") <--> MAX(" << Xmax << "," << Ymax << "," << Zmax 
			      << ")" << std::endl;*/

	// --------------------------------------------------------------------
	//			STORING EDGES
	// --------------------------------------------------------------------

	gridfile >> line;		// it should be: '}'
	gridfile >> line; 		// it should be: 'Edges'
	gridfile >> line; 		// it should be: '(nb_edges)'
	line.erase(line.begin());
	line.erase(line.end()-1);	
	nb = atoi(line.c_str());	// it should be now 'nb_edges'


	if (nb != nb_edges) { 
		ERRORMSG("Number of edges mismatch !! : " << nb << " != " << nb_edges);
		/*std::cout << "ERROR - Nb of edges mismatch!!: " << nb
			  << " != " << nb_edges << std::endl;*/
	};

	gridfile >> line; // it should be ')'

	for (unsigned i=0; i<nb_edges; i++){
		gridfile >> line;	Edges[i].first = atoi(line.c_str());
		gridfile >> line;	Edges[i].second = atoi(line.c_str());
	};
	LOWMSG(nb_edges << "\t edges stored.");
	//std::cout << nb_edges << "\t edges stored." << std::endl;

	// --------------------------------------------------------------------
	//			STORING FACES
	// --------------------------------------------------------------------

	gridfile >> line; 	// it should be: '}'
	gridfile >> line; 	// it should be: 'Faces'
	gridfile >> line; 	// it should be: '(nb_faces)'
	line.erase(line.begin());
	line.erase(line.end()-1);
	nb = atoi(line.c_str());	// it should be now 'nb_faces'


	if (nb != nb_faces) { 
		ERRORMSG("Number of faces mismatch !! : " << nb << " != " << nb_faces);
		/*std::cout << "ERRORMSG(" - Nb of faces mismatch!!: " << nb
			  << " != " << nb_faces << std::endl;*/
	};

	gridfile >> line; 	// it should be ')'


	for (unsigned i=0; i<nb_faces; i++){
		gridfile >> face_code; Faces[i].code = face_code;
		// ALERT!!! : mind the sign later, as it's also stored!!
		switch (face_code) {
			case 3:		/* TRIANGLE */ {
				gridfile >> line;	Faces[i].edge[0] = atoi(line.c_str());
				gridfile >> line;	Faces[i].edge[1] = atoi(line.c_str());
				gridfile >> line;	Faces[i].edge[2] = atoi(line.c_str());
				break; }
			case 4:		/* SQUARE */ {
				gridfile >> line;	Faces[i].edge[0] = atoi(line.c_str());
				gridfile >> line;	Faces[i].edge[1] = atoi(line.c_str());
				gridfile >> line;	Faces[i].edge[2] = atoi(line.c_str());
				gridfile >> line;	Faces[i].edge[3] = atoi(line.c_str());
				break; }
			default: 
				WARNINGMSG("Code of Face(" << i << ") = "<< face_code <<" not recognized !!");
				//std::cout << "Code of Face(" << i << ") = "<< face_code <<" not recognized !!" << std::endl;
		}
	};
	LOWMSG(nb_faces << "\t faces stored.");
	//std::cout << nb_faces << "\t faces stored." << std::endl;

	// --------------------------------------------------------------------
	//			STORING ELEMENTS
	// --------------------------------------------------------------------

	gridfile >> line; 	// it should be: '}'

	while (line != "Elements"){
		gridfile >> line;
	};

	gridfile >> line; 	// it should be: '(nb_elements)' << std::endl
	line.erase(line.begin());
	line.erase(line.end()-1);
	nb = atoi(line.c_str());	// it should be now 'nb_elements'


	if (nb != nb_elements) { 
		ERRORMSG("Number of elements mismatch !! : " << nb << " != " << nb_elements);
	};

	gridfile >> line; 	// it should be ')'
	//gridfile >> line; 	// it should be: '{'
	nb_tetrahedra = nb_cubes = nb_triangles = 0;

	for (unsigned i=0; i<nb_elements; i++){
		gridfile >> elem_code; Elements[i].code = elem_code;
		switch (elem_code) {
			case 2:		/* TRIANGLE */ {
				gridfile >> line; 	Elements[i].face[0] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[1] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[2] = atoi(line.c_str());
				nb_triangles++; break; }
			case 5:		/* TETRAHEDRON */ {
				gridfile >> line;	Elements[i].face[0] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[1] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[2] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[3] = atoi(line.c_str());
				nb_tetrahedra++; break; }
			case 8:		/* CUBE */ {
				gridfile >> line;	Elements[i].face[0] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[1] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[2] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[3] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[4] = atoi(line.c_str());
				gridfile >> line;	Elements[i].face[5] = atoi(line.c_str());
				nb_cubes++; break; }
			default:	
				WARNINGMSG("Code of Element(" << i << ") = "<< elem_code <<" not recognized !!");
				//std::cout << "Code of Element(" << i << ") = "<< elem_code <<" not recognized !!" << std::endl;
		}


	};
	LOWMSG(nb_elements << "\t elements stored.");
	//std::cout << nb_elements << "\t elements stored." << std::endl;


	// --------------------------------------------------------------------
	//			STORING MATERIALS
	// --------------------------------------------------------------------

	gridfile >> line;	// } //
	gridfile >> end;

	while (end != "}") {
		// else: end = Region
		gridfile >> line; gridfile >> line; gridfile >> line; gridfile >> line;
		// ("RegionName") // { // material // = //
		gridfile >> Material;
		gridfile >> line; gridfile >> line; gridfile >> line;
		//Elements // (NumberOfElementsInThatRegion) // { //
		gridfile >> line;
		while (line != "}") {
			Elements[atoi(line.c_str())].material = Material;
			gridfile >> line;
		}
		gridfile >> line; // } //
		gridfile >> end; // } ? Region
	};
	LOWMSG("MeshParser file successfully processed.\nPreparing uniform grid from valuable information ...");
	//std::cout << "MeshParser file successfully processed. Calculating valuable information." << std::endl;
	gridfile.close();
	
	//std::cout << "Gridfile closed " << std::endl;

	// --------------------------------------------------------------------
	//		  STORING VERTICES INTO ELEMENTS
	// --------------------------------------------------------------------



	//logfile << " Starting lopp through all the faces (nb_faces = " << nb_faces << " ) " << endl;
	for (unsigned i=0; i<nb_faces; i++)
	{
		face_code = Faces[i].code;
		switch (face_code) {
			case 3: 	/* TRIANGLE */ {
				edge_id = Faces[i].edge[0];
				if (edge_id < 0) { edge_id = (-1)*edge_id - 1;};

				Faces[i].vertex[0] = Edges[edge_id].first;
				Faces[i].vertex[1] = Edges[edge_id].second;


				edge_id = Faces[i].edge[1];
				if (edge_id < 0) { edge_id = (-1)*edge_id - 1;};

				if (Edges[edge_id].first == Faces[i].vertex[0] ||
				    Edges[edge_id].first == Faces[i].vertex[1] ) { 
					Faces[i].vertex[2] = Edges[edge_id].second;
				} else {
				 	Faces[i].vertex[2] = Edges[edge_id].first;
				}; break; }
			case 4: 	/* SQUARE */ {
				edge_id = Faces[i].edge[0];
				if (edge_id < 0) { edge_id = (-1)*edge_id - 1;};

				Faces[i].vertex[0] = Edges[edge_id].first;
				Faces[i].vertex[1] = Edges[edge_id].second;


				edge_id = Faces[i].edge[1];
				if (edge_id < 0) { edge_id = (-1)*edge_id - 1;};

				if (Edges[edge_id].first == Faces[i].vertex[0] ||
					Edges[edge_id].first == Faces[i].vertex[1] ) { 
					Faces[i].vertex[2] = Edges[edge_id].second;
				} else {
				 	Faces[i].vertex[2] = Edges[edge_id].first;
				};


				edge_id = Faces[i].edge[2];
				if (edge_id < 0) { edge_id = (-1)*edge_id - 1;};
				if (Edges[edge_id].first == Faces[i].vertex[0] ||
					Edges[edge_id].first == Faces[i].vertex[1] ||
					Edges[edge_id].first == Faces[i].vertex[2]) { 
					Faces[i].vertex[3] = Edges[edge_id].second;
				} else {
				 	Faces[i].vertex[3] = Edges[edge_id].first;
				}; break; }
			default:
				ERRORMSG("Code of Facet(" << i
					<< ") = "<< face_code <<" not recognized !!");
		}

		/*logfile << " Face [" << i << "] has vertices : ( " << Faces[i].vertices[0] <<
			" , " << Faces[i].vertices[1] << " , " << Faces[i].vertices[2] << " ) " << endl;*/

	}; // Close loop for all Faces

	//std::cout << "Faces are now also defined using vertices." << std::endl;


	LOWMSG("Starting loop through "	<< nb_elements << " elements:"
			<< "\n\t\tTetrahedra\t= "	<< nb_tetrahedra
			<< "\n\t\tCubes\t\t= "	<< nb_cubes
			<< "\n\t\tTriangles\t= "	<< nb_triangles);
	/*std::cout	<< "Starting loop through "	<< nb_elements << " elements:"
				<< "\n\t\tTetrahedra\t= "	<< nb_tetrahedra
				<< "\n\t\tCubes\t\t= "	<< nb_cubes
				<< "\n\t\tTriangles\t= "	<< nb_triangles << std::endl;*/


	for (unsigned i=0; i<nb_elements; i++)
	{
		elem_code = Elements[i].code;
		switch (elem_code) {
			case 2:		/* TRIANGLE */ {
				WARNINGMSG("2D Element detected.");
				break; }
			case 5:		/* TETRAHEDRON */ {
				face_id = Elements[i].face[0];
				if (face_id < 0) { face_id = (-1)*face_id -1 ; };

				Elements[i].vertex[0] = Faces[face_id].vertex[0];
				Elements[i].vertex[1] = Faces[face_id].vertex[1];
				Elements[i].vertex[2] = Faces[face_id].vertex[2];

				face_id = Elements[i].face[1]; // it might be also 2 or 3
				if (face_id < 0) { face_id = (-1)*face_id - 1; };

				if (  Faces[face_id].vertex[0] != Elements[i].vertex[0] &&
					  Faces[face_id].vertex[0] != Elements[i].vertex[1] &&
					  Faces[face_id].vertex[0] != Elements[i].vertex[2])
				      { Elements[i].vertex[3] = Faces[face_id].vertex[0]; }
				else if ( Faces[face_id].vertex[1] != Elements[i].vertex[0] &&
					  Faces[face_id].vertex[1] != Elements[i].vertex[1] &&
					  Faces[face_id].vertex[1] != Elements[i].vertex[2])
				      { Elements[i].vertex[3] = Faces[face_id].vertex[1]; }
				else  { Elements[i].vertex[3] = Faces[face_id].vertex[2]; };
				break; }
			case 8:		/* CUBE */ {
				int vx = 3;
				face_id = Elements[i].face[0];
				if (face_id < 0) { face_id = (-1)*face_id -1 ; };

				Elements[i].vertex[0] = Faces[face_id].vertex[0];
				Elements[i].vertex[1] = Faces[face_id].vertex[1];
				Elements[i].vertex[2] = Faces[face_id].vertex[2];
				Elements[i].vertex[3] = Faces[face_id].vertex[3];


				face_id = Elements[i].face[1];
				if (face_id < 0) { face_id = (-1)*face_id - 1; };

				if(  Faces[face_id].vertex[0] != Elements[i].vertex[0] &&
				     Faces[face_id].vertex[0] != Elements[i].vertex[1] &&
				     Faces[face_id].vertex[0] != Elements[i].vertex[2] &&
				     Faces[face_id].vertex[0] != Elements[i].vertex[3])
				   { vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[0]; }
				if ( Faces[face_id].vertex[1] != Elements[i].vertex[0] &&
				     Faces[face_id].vertex[1] != Elements[i].vertex[1] &&
				     Faces[face_id].vertex[1] != Elements[i].vertex[2] &&
				     Faces[face_id].vertex[1] != Elements[i].vertex[3])
				   { vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[1]; }
				if ( Faces[face_id].vertex[2] != Elements[i].vertex[0] &&
				     Faces[face_id].vertex[2] != Elements[i].vertex[1] &&
				     Faces[face_id].vertex[2] != Elements[i].vertex[2] &&
				     Faces[face_id].vertex[2] != Elements[i].vertex[3])
				   { vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[2]; }
				if(  Faces[face_id].vertex[3] != Elements[i].vertex[0] &&
				     Faces[face_id].vertex[3] != Elements[i].vertex[1] &&
				     Faces[face_id].vertex[3] != Elements[i].vertex[2] &&
				     Faces[face_id].vertex[3] != Elements[i].vertex[3])
				   { vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[3]; }

				if (vx == 5) {
					face_id = Elements[i].face[2];
					if (face_id < 0) { face_id = (-1)*face_id - 1; };
					if(	  Faces[face_id].vertex[0] != Elements[i].vertex[0] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[1] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[2] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[3] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[4] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[5])
					      {	vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[0]; }
					else if ( Faces[face_id].vertex[1] != Elements[i].vertex[0] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[1] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[2] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[3] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[4] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[5])
					      {	vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[1]; }
					else if ( Faces[face_id].vertex[2] != Elements[i].vertex[0] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[1] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[2] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[3] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[4] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[5])
					      {	vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[2]; }
					else if ( Faces[face_id].vertex[3] != Elements[i].vertex[0] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[1] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[2] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[3] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[4] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[5])
					      {	vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[3]; };
				};
				if (vx == 6) {
					face_id = Elements[i].face[3];
					if (face_id < 0) { face_id = (-1)*face_id - 1; };
					if(	  Faces[face_id].vertex[0] != Elements[i].vertex[0] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[1] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[2] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[3] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[4] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[5] &&
						  Faces[face_id].vertex[0] != Elements[i].vertex[6])
					      {	vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[0]; }
					else if ( Faces[face_id].vertex[1] != Elements[i].vertex[0] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[1] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[2] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[3] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[4] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[5] &&
						  Faces[face_id].vertex[1] != Elements[i].vertex[6])
					      {	vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[1]; }
					else if ( Faces[face_id].vertex[2] != Elements[i].vertex[0] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[1] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[2] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[3] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[4] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[5] &&
						  Faces[face_id].vertex[2] != Elements[i].vertex[6])
					      {	vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[2]; }
					else if ( Faces[face_id].vertex[3] != Elements[i].vertex[0] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[1] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[2] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[3] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[4] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[5] &&
						  Faces[face_id].vertex[3] != Elements[i].vertex[6])
					      {	vx++; Elements[i].vertex[vx] = Faces[face_id].vertex[3]; }
				}; break; }
			default:
				ERRORMSG("Code of Element(" << i
					<< ") = "<< elem_code <<" not recognized !!");
		}
	}; // Close loop for all Elements

	//std::cout << "Elements are now also defined using vertices." << std::endl;

	// --------------------------------------------------------------------
	//		  	CREATION OF OCTREE
	// --------------------------------------------------------------------


	// Resizing Octree
	Octree.resize(HULL); // X-axis
	for (unsigned i=0; i<HULL; i++) {
		Octree[i].resize(HULL); // Y-axis
		for (unsigned j=0; j<HULL; j++) {
			Octree[i][j].resize(HULL); // Z-axis
		}
	}



	// Loop along all the Elements

	for (unsigned elem_id=0; elem_id<nb_elements; elem_id++) {
		elem_code = Elements[elem_id].code;

		switch (elem_code) {
			case 2:		/* TRIANGLE */ {
				WARNINGMSG("2D Element detected.");
				break; }
			case 5:		/* TETRAHEDRON */ {

				vertex_id = Elements[elem_id].vertex[0]; // Vertex #1  
				xmin = xmax = Vertices[vertex_id].x;
				ymin = ymax = Vertices[vertex_id].y;
				zmin = zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[1]; // Vertex #2
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[2]; // Vertex #3
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[3]; // Vertex #4
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;
				break; }
			case 8:		/* CUBE */ {

				vertex_id = Elements[elem_id].vertex[0]; // Vertex #1  
				xmin = xmax = Vertices[vertex_id].x;
				ymin = ymax = Vertices[vertex_id].y;
				zmin = zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[1]; // Vertex #2
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[2]; // Vertex #3
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[3]; // Vertex #4
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[4]; // Vertex #5
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[5]; // Vertex #6
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[6]; // Vertex #7
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;

				vertex_id = Elements[elem_id].vertex[7]; // Vertex #8
				if (Vertices[vertex_id].x < xmin) xmin = Vertices[vertex_id].x;
				if (Vertices[vertex_id].x > xmax) xmax = Vertices[vertex_id].x;
				if (Vertices[vertex_id].y < ymin) ymin = Vertices[vertex_id].y;
				if (Vertices[vertex_id].y > ymax) ymax = Vertices[vertex_id].y;
				if (Vertices[vertex_id].z < zmin) zmin = Vertices[vertex_id].z;
				if (Vertices[vertex_id].z > zmax) zmax = Vertices[vertex_id].z;
				break; }
			default:
				ERRORMSG("Code of Element(" << elem_id
					<< ") = "<< elem_code <<" not recognized !!");
		} // close switch
				

		// Coordinates of the diagonal vertices of the isothetic box hull

		Hull[0].x = xmin;	Hull[1].x = xmax;
		Hull[0].y = ymin;	Hull[1].y = ymax;
		Hull[0].z = zmin;	Hull[1].z = zmax;


		// Find intersecting cubes with box hull diagonal vertices
		for (int h=0; h<2; h++) { 
			x = (Hull[h].x + Xshift) / (Xmax - Xmin) * HULL;
			y = (Hull[h].y + Yshift) / (Ymax - Ymin) * HULL;
			z = (Hull[h].z + Zshift) / (Zmax - Zmin) * HULL;
			
			if (x == HULL ) x = HULL - 1;
			if (y == HULL ) y = HULL - 1;
			if (z == HULL ) z = HULL - 1;
			idx = (int) x;
			idy = (int) y;
			idz = (int) z;
			
			Cube[h].first  = idx;
			Cube[h].second = idy;
			Cube[h].third  = idz;
			
		};
		// Are they the same cube? Yes -> Store Element in Octree[idx][idy][idz]
		if ((Cube[0].first == Cube[1].first) && (Cube[0].second == Cube[1].second) && (Cube[0].third == Cube[1].third) ) { 
			Octree.at(Cube[0].first).at(Cube[0].second).at(Cube[0].third).push_back(elem_id); 
		}else {
		// Find cubes within the volume of cubes at hull vertices and attach them elem_id
			for (int ix=Cube[0].first; ix<=Cube[1].first; ix++){
				for (int iy=Cube[0].second; iy<=Cube[1].second; iy++){
					for (int iz=Cube[0].third; iz<=Cube[1].third; iz++){
						Octree.at(ix).at(iy).at(iz).push_back(elem_id); 
					}; // end FOR
				}; // end FOR
			}; // end FOR
	
		}; // end ELSE
	}; // End of nb_elements loop

	LOWMSG("Uniform grid successfully created and stored. Grid information:");
	//std::cout << "Uniform grid successfully created and stored. Grid information:" << std::endl;

	// Find out max number of tetra_id in hulls
	unsigned maxx = 0;
	unsigned minn = 1E6;
	unsigned objects = 0;
	for (unsigned i=0; i<HULL; i++) {
		for (unsigned j=0; j<HULL; j++) {
			for (unsigned k=0; k<HULL; k++) {
				if (Octree[i][j][k].size() > maxx) maxx = Octree[i][j][k].size();
				if (Octree[i][j][k].size() < minn) minn = Octree[i][j][k].size();
				objects += Octree[i][j][k].size();
			}
		}
	}

	//std::cout << "Vertices size: " << (float)sizeof(Vertices)/1024.0/1024.0  << " MB" << std::endl;
	//std::cout << "Elements size: " << (float)sizeof(Elements)/1024.0/1024.0  << " MB" << std::endl;
	/*std::cout << "\tSize in memory: " << (float)sizeof(int)*objects/1024.0/1024.0 << " MB" << std::endl;
	std::cout << "\tContent information:" << std::endl;
	std::cout << "\t\tMaximum size\t= " 	<< maxx    << std::endl;
	std::cout << "\t\tMinimum size\t= " 	<< minn    << std::endl;
	std::cout << "\t\tTotal elements\t= "	<< objects << std::endl;*/
	LOWMSG("\tSize in memory: " << (float)sizeof(int)*objects/1024.0/1024.0 << " MB");
	LOWMSG("\tContent information:");
	LOWMSG("\t\tMaximum size\t= " 	<< maxx);
	LOWMSG("\t\tMinimum size\t= " 	<< minn);
	LOWMSG("\t\tTotal elements\t= "	<< objects);

}

// -----------------------------------------------------------------------------
// Find Octree containing p(x,y,z)
// -----------------------------------------------------------------------------

MeshParser::int3 MeshParser::FindOctree(float x, float y, float z) const
{
	int3 ret;
	ret.first  = (int)( (x + Xshift) / dimX * HULL );
	ret.second = (int)( (y + Yshift) / dimY * HULL );
	ret.third  = (int)( (z + Zshift) / dimZ * HULL );
	if (ret.first  == HULL ) ret.first  = HULL - 1;
	if (ret.second == HULL ) ret.second = HULL - 1;
	if (ret.third  == HULL ) ret.third  = HULL - 1;

	return ret;
}

// -----------------------------------------------------------------------------
// Calculate the barycentric coordinates of p(x,y,z) respect to tetraID
// -----------------------------------------------------------------------------

MeshParser::float3 MeshParser::GetBarycentricCoords(float x, float y, float z, int tetraID) const
{

	// Determinant det = Q4 (Q2 x Q3) where Q2 = Vertex2-Vertex1, Q3 = Vertex3-Vertex1, Q4 = Vertex4-Vertex1.
	//
	float det;
	det = 	 (Vertices[Elements[tetraID].vertex[3]].x - Vertices[Elements[tetraID].vertex[0]].x)*
		((Vertices[Elements[tetraID].vertex[1]].y - Vertices[Elements[tetraID].vertex[0]].y)*
		 (Vertices[Elements[tetraID].vertex[2]].z - Vertices[Elements[tetraID].vertex[0]].z)-
		 (Vertices[Elements[tetraID].vertex[1]].z - Vertices[Elements[tetraID].vertex[0]].z)*
		 (Vertices[Elements[tetraID].vertex[2]].y - Vertices[Elements[tetraID].vertex[0]].y))
	       + (Vertices[Elements[tetraID].vertex[3]].y - Vertices[Elements[tetraID].vertex[0]].y)*
		((Vertices[Elements[tetraID].vertex[1]].z - Vertices[Elements[tetraID].vertex[0]].z)*
		 (Vertices[Elements[tetraID].vertex[2]].x - Vertices[Elements[tetraID].vertex[0]].x)-
		 (Vertices[Elements[tetraID].vertex[1]].x - Vertices[Elements[tetraID].vertex[0]].x)*
		 (Vertices[Elements[tetraID].vertex[2]].z - Vertices[Elements[tetraID].vertex[0]].z)) 
	       + (Vertices[Elements[tetraID].vertex[3]].z - Vertices[Elements[tetraID].vertex[0]].z)*
		((Vertices[Elements[tetraID].vertex[1]].x - Vertices[Elements[tetraID].vertex[0]].x)*
		 (Vertices[Elements[tetraID].vertex[2]].y - Vertices[Elements[tetraID].vertex[0]].y)-
		 (Vertices[Elements[tetraID].vertex[1]].y - Vertices[Elements[tetraID].vertex[0]].y)*
		 (Vertices[Elements[tetraID].vertex[2]].x - Vertices[Elements[tetraID].vertex[0]].x));

	float3 ret;
	ret.x = ( (x - Vertices[Elements[tetraID].vertex[0]].x)*
			((Vertices[Elements[tetraID].vertex[2]].y - Vertices[Elements[tetraID].vertex[0]].y)*
			 (Vertices[Elements[tetraID].vertex[3]].z - Vertices[Elements[tetraID].vertex[0]].z)-
			 (Vertices[Elements[tetraID].vertex[2]].z - Vertices[Elements[tetraID].vertex[0]].z)*
			 (Vertices[Elements[tetraID].vertex[3]].y - Vertices[Elements[tetraID].vertex[0]].y)) 
	     			+ (y - Vertices[Elements[tetraID].vertex[0]].y)*
			((Vertices[Elements[tetraID].vertex[2]].z - Vertices[Elements[tetraID].vertex[0]].z)*
			 (Vertices[Elements[tetraID].vertex[3]].x - Vertices[Elements[tetraID].vertex[0]].x)-
			 (Vertices[Elements[tetraID].vertex[2]].x - Vertices[Elements[tetraID].vertex[0]].x)*
			 (Vertices[Elements[tetraID].vertex[3]].z - Vertices[Elements[tetraID].vertex[0]].z)) 
	     			+ (z - Vertices[Elements[tetraID].vertex[0]].z)*
			((Vertices[Elements[tetraID].vertex[2]].x - Vertices[Elements[tetraID].vertex[0]].x)*
			 (Vertices[Elements[tetraID].vertex[3]].y - Vertices[Elements[tetraID].vertex[0]].y)-
			 (Vertices[Elements[tetraID].vertex[2]].y - Vertices[Elements[tetraID].vertex[0]].y)*
			 (Vertices[Elements[tetraID].vertex[3]].x - Vertices[Elements[tetraID].vertex[0]].x)) )/det;
	ret.y = ( (x - Vertices[Elements[tetraID].vertex[0]].x)*
			((Vertices[Elements[tetraID].vertex[3]].y - Vertices[Elements[tetraID].vertex[0]].y)*
			 (Vertices[Elements[tetraID].vertex[1]].z - Vertices[Elements[tetraID].vertex[0]].z)-
			 (Vertices[Elements[tetraID].vertex[3]].z - Vertices[Elements[tetraID].vertex[0]].z)*
			 (Vertices[Elements[tetraID].vertex[1]].y - Vertices[Elements[tetraID].vertex[0]].y)) 
	     			+ (y - Vertices[Elements[tetraID].vertex[0]].y)*
			((Vertices[Elements[tetraID].vertex[3]].z - Vertices[Elements[tetraID].vertex[0]].z)*
			 (Vertices[Elements[tetraID].vertex[1]].x - Vertices[Elements[tetraID].vertex[0]].x)-
			 (Vertices[Elements[tetraID].vertex[3]].x - Vertices[Elements[tetraID].vertex[0]].x)*
			 (Vertices[Elements[tetraID].vertex[1]].z - Vertices[Elements[tetraID].vertex[0]].z)) 
	     			+ (z - Vertices[Elements[tetraID].vertex[0]].z)*
			((Vertices[Elements[tetraID].vertex[3]].x - Vertices[Elements[tetraID].vertex[0]].x)*
			 (Vertices[Elements[tetraID].vertex[1]].y - Vertices[Elements[tetraID].vertex[0]].y)-
			 (Vertices[Elements[tetraID].vertex[3]].y - Vertices[Elements[tetraID].vertex[0]].y)*
			 (Vertices[Elements[tetraID].vertex[1]].x - Vertices[Elements[tetraID].vertex[0]].x)) )/det;
	ret.z = ( (x - Vertices[Elements[tetraID].vertex[0]].x)*
			((Vertices[Elements[tetraID].vertex[1]].y - Vertices[Elements[tetraID].vertex[0]].y)*
			 (Vertices[Elements[tetraID].vertex[2]].z - Vertices[Elements[tetraID].vertex[0]].z)-
			 (Vertices[Elements[tetraID].vertex[1]].z - Vertices[Elements[tetraID].vertex[0]].z)*
			 (Vertices[Elements[tetraID].vertex[2]].y - Vertices[Elements[tetraID].vertex[0]].y)) 
	     			+ (y - Vertices[Elements[tetraID].vertex[0]].y)*
			((Vertices[Elements[tetraID].vertex[1]].z - Vertices[Elements[tetraID].vertex[0]].z)*
			 (Vertices[Elements[tetraID].vertex[2]].x - Vertices[Elements[tetraID].vertex[0]].x)-
			 (Vertices[Elements[tetraID].vertex[1]].x - Vertices[Elements[tetraID].vertex[0]].x)*
			 (Vertices[Elements[tetraID].vertex[2]].z - Vertices[Elements[tetraID].vertex[0]].z)) 
	     			+ (z - Vertices[Elements[tetraID].vertex[0]].z)*
			((Vertices[Elements[tetraID].vertex[1]].x - Vertices[Elements[tetraID].vertex[0]].x)*
			 (Vertices[Elements[tetraID].vertex[2]].y - Vertices[Elements[tetraID].vertex[0]].y)-
			 (Vertices[Elements[tetraID].vertex[1]].y - Vertices[Elements[tetraID].vertex[0]].y)*
			 (Vertices[Elements[tetraID].vertex[2]].x - Vertices[Elements[tetraID].vertex[0]].x)) )/det;


	return ret;
}


// -----------------------------------------------------------------------------
// Find element (inside Octree) containing p(x,y,z)
// -----------------------------------------------------------------------------

int MeshParser::FindElement(float x, float y, float z, int3 oct) const // Registers: 3x int + 1x float3 + 2x printf
{
	int Eid,tempid;
	tempid = Eid = -1;
	for (unsigned i=0; i<Octree[oct.first][oct.second][oct.third].size(); i++){
		Eid = Octree[oct.first][oct.second][oct.third][i];
		switch (Elements[Eid].code) {
			case (2): /* TRIANGLE */ {
				WARNINGMSG("Material requested from 2D element.");
				break; }
			case (5): /* TETRAHEDRON */ {
				float3 bary;
				bary = GetBarycentricCoords(x,y,z,Eid);
				if(bary.x >= 0 && bary.y >= 0 && bary.z >= 0){
					if(bary.x+bary.y+bary.z <= 1.000001){
						return Eid;
					}else if (bary.x+bary.y+bary.z > 1.000001 && bary.x+bary.y+bary.z <= 1.05){
						tempid = Eid;
					}
				} 
				break; }
			case (8): /* CUBE */ {
				int2 diagonal;
				diagonal = getDiagonal(Eid);
				if (Vertices[diagonal.first].x <= x && Vertices[diagonal.second].x >= x &&
				    Vertices[diagonal.first].y <= y && Vertices[diagonal.second].y >= y &&
				    Vertices[diagonal.first].z <= z && Vertices[diagonal.second].z >= z)
					return Eid;

				break; }
			default:
				ERRORMSG("Code of Element(" << i
					  << ") = "<< elem_code <<" not recognized !!");
			}

	}	

	if(tempid >= 0){ 
		WARNINGMSG("Octree[" << oct.first <<"][" << oct.second << "][" << oct.third << "] barely contains (x,y,z) = "
							 << x << " , " << y << " , " << z);
		return tempid;
	}
	ERRORMSG("Any Element from Octree[" <<oct.first << "][" << oct.second << "][" << oct.third <<"] contains (x,y,z) = "
			<< x <<" , "<< y <<" , "<< z);
	return Eid;
}

// -----------------------------------------------------------------------------
// Get diagonal (== ID of two vertices) of cubic element EID
// -----------------------------------------------------------------------------

MeshParser::int2 MeshParser::getDiagonal(int EID) const
{
	
	int2 point;
	point.first  = Elements[EID].vertex[0];
	point.second = Elements[EID].vertex[0];
	for (int n=1; n<8; n++) {
		if( Vertices[Elements[EID].vertex[n]].x <= Vertices[point.first].x &&
		    Vertices[Elements[EID].vertex[n]].y <= Vertices[point.first].y &&
		    Vertices[Elements[EID].vertex[n]].z <= Vertices[point.first].z ) 
			point.first = Elements[EID].vertex[n];
		if( Vertices[Elements[EID].vertex[n]].x >= Vertices[point.second].x &&
		    Vertices[Elements[EID].vertex[n]].y >= Vertices[point.second].y &&
		    Vertices[Elements[EID].vertex[n]].z >= Vertices[point.second].z ) 
			point.second = Elements[EID].vertex[n];
	}
	return point;
}

std::string MeshParser::getMaterial(float X, float Y, float Z) const
{
	int3 OCT = FindOctree(X,Y,Z);
	int elem = FindElement(X,Y,Z,OCT);
	return Elements[elem].material;
}



