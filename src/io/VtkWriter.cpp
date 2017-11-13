/**
* \file VtkWriter.cpp
* \brief Vtk writer interface.
* \date December 01, 2014
* \author Benoit SKlenard
* \bug No known bugs.
* 
* - Write binary files not implemented
* - Only ImageData available (VtkWriterImageData class)
*/

#include <cassert>
#include <fstream>
#include <ios>
#include <iomanip>
#include <iostream>

#include "Base64.h"
#include "VtkWriter.h"

#define VTK_VERTEX         1
#define VTK_LINE           3 
#define VTK_POLY_VERTEX    2 
#define VTK_POLY_LINE      4 
#define VTK_TRIANGLE_STRIP 6 
#define VTK_TRIANGLE       5 
#define VTK_QUAD           9 
#define VTK_POLYGON        7 
#define VTK_PIXEL          8 
#define VTK_TETRA          10
#define VTK_VOXEL          11
#define VTK_HEXAHEDRON     12
 
namespace IO {

const std::string VtkWriter::vtkTypeString[] = { "Int8", "UInt8", "Int16", "UInt16", "Int32", "UInt32", "Int64", "UInt64", "Float32", "Float64" };
const std::string VtkWriter::vtkFormatString[] = { "ascii", "binary" };

VtkWriter::VtkWriter() {
    // test byte order
    short int number = 0x1;
    char *numPtr = (char*)&number;

    _endianness = (numPtr[0] == 1) ? "LittleEndian" : "BigEndian";
}

void VtkWriter::indent(std::ostream &os, int shift) { 
    for (int i=0; i < shift; ++i)
        os << " ";
}

template <typename T>
void VtkWriter::writeBinaryScalars(std::ostream &os, const std::vector<T> &v) { 
    if (v.size()) {
		uint32_t len = v.size()*sizeof(T);

        std::string b64_len;
        std::string b64_data;

        base64_encode(reinterpret_cast<const char *>(&len), sizeof(uint32_t), b64_len);
        base64_encode(reinterpret_cast<const char *>(&v[0]), len, b64_data);

		os.write(b64_len.c_str(), b64_len.size());
		os.write(b64_data.c_str(), b64_data.size());
    }
}


template <typename T>
void VtkWriter::writeAsciiScalars(std::ostream &os, const std::vector<T> &v, int shift) {
    int n = v.size();
    if (n) {
        indent(os,shift);   
        for (int i=0; i<n; ++i)
            os << v[i] << " ";
        os << "\n";
    }
}

template <typename T>
void VtkWriter::writeAsciiVectors(std::ostream &os, const std::vector<std::vector<T> > &v, int nComponents, int shift) {
    int n = v.size();
    for (int i=0; i<n; ++i) {
        indent(os,shift);
        for (int j=0; j<nComponents; ++j) {
            os << v[i][j] << " ";
        }
        os << "\n";
    } 
}

void VtkWriter::writeScalarDataArrayBlock(std::ostream &os, const DataSet &set, int shift)
{
    indent(os, shift);
    os << "<DataArray Name=\"" << set._name << "\" " \
                     "type=\"" << vtkTypeString[set._type] << "\" "\
                     "format=\"" << vtkFormatString[set._format] << "\">\n";
    if (set._format == ASCII) {
    	switch (set._type) {
        	case INT8    : writeAsciiScalars(os, *((std::vector<int8_t> *)set._data), shift+2); break;
        	case UINT8   : writeAsciiScalars(os, *((std::vector<uint8_t> *)set._data), shift+2); break;
        	case INT16   : writeAsciiScalars(os, *((std::vector<int16_t> *)set._data), shift+2); break;
        	case UINT16  : writeAsciiScalars(os, *((std::vector<uint16_t> *)set._data), shift+2); break;
	        case INT32   : writeAsciiScalars(os, *((std::vector<int32_t> *)set._data), shift+2); break;
        	case UINT32  : writeAsciiScalars(os, *((std::vector<uint32_t> *)set._data), shift+2); break;
    	    case INT64   : writeAsciiScalars(os, *((std::vector<int64_t> *)set._data), shift+2); break;
	        case UINT64  : writeAsciiScalars(os, *((std::vector<uint64_t> *)set._data), shift+2); break;
        	case FLOAT32 : writeAsciiScalars(os, *((std::vector<float> *)set._data), shift+2); break;
    	    case FLOAT64 : writeAsciiScalars(os, *((std::vector<double> *)set._data), shift+2); break;
			indent(os, shift);
	    }
	}
	else {
    	switch (set._type) {
        	case INT8    : writeBinaryScalars(os, *((std::vector<int8_t> *)set._data)); break;
        	case UINT8   : writeBinaryScalars(os, *((std::vector<uint8_t> *)set._data)); break;
        	case INT16   : writeBinaryScalars(os, *((std::vector<int16_t> *)set._data)); break;
        	case UINT16  : writeBinaryScalars(os, *((std::vector<uint16_t> *)set._data)); break;
	        case INT32   : writeBinaryScalars(os, *((std::vector<int32_t> *)set._data)); break;
        	case UINT32  : writeBinaryScalars(os, *((std::vector<uint32_t> *)set._data)); break;
    	    case INT64   : writeBinaryScalars(os, *((std::vector<int64_t> *)set._data)); break;
	        case UINT64  : writeBinaryScalars(os, *((std::vector<uint64_t> *)set._data)); break;
        	case FLOAT32 : writeBinaryScalars(os, *((std::vector<float> *)set._data)); break;
    	    case FLOAT64 : writeBinaryScalars(os, *((std::vector<double> *)set._data)); break;
	    }
	}
    os << "</DataArray>\n";
}

void VtkWriter::writeVectorDataArrayBlock(std::ostream &os, const DataSet &set, int shift)
{
    indent(os, shift);
    os << "<DataArray Name=\"" << set._name << "\" " \
                     "NumberOfComponents=\"" << set._nComponents << "\" " \
                     "type=\"" << vtkTypeString[set._type] << "\" " \
                     "format=\"" << vtkFormatString[set._format] << "\">\n";
    switch (set._type) {
        case INT8    : writeAsciiVectors(os, *((std::vector<std::vector<int8_t> > *)set._data), set._nComponents, shift+2); break;
        case UINT8   : writeAsciiVectors(os, *((std::vector<std::vector<uint8_t> > *)set._data), set._nComponents, shift+2); break;
        case INT16   : writeAsciiVectors(os, *((std::vector<std::vector<int16_t> > *)set._data), set._nComponents, shift+2); break;
        case UINT16  : writeAsciiVectors(os, *((std::vector<std::vector<uint16_t> > *)set._data), set._nComponents, shift+2); break;
        case INT32   : writeAsciiVectors(os, *((std::vector<std::vector<int32_t> > *)set._data), set._nComponents, shift+2); break;
        case UINT32  : writeAsciiVectors(os, *((std::vector<std::vector<uint32_t> > *)set._data), set._nComponents, shift+2); break;
        case INT64   : writeAsciiVectors(os, *((std::vector<std::vector<int64_t> > *)set._data), set._nComponents, shift+2); break;
        case UINT64  : writeAsciiVectors(os, *((std::vector<std::vector<uint64_t> > *)set._data), set._nComponents, shift+2); break;
        case FLOAT32 : writeAsciiVectors(os, *((std::vector<std::vector<float> > *)set._data), set._nComponents, shift+2); break;
        case FLOAT64 : writeAsciiVectors(os, *((std::vector<std::vector<double> > *)set._data), set._nComponents, shift+2); break;
    }
    indent(os, shift);
    os << "</DataArray>\n";
}

void VtkWriter::write() {
    std::ofstream ofs;

    ofs.open(_filename.c_str(), std::ofstream::out | std::ofstream::binary);
    write(ofs);
    ofs.close();
}

void VtkWriter::writePointDataBlock(std::ostream &os)
{
    os << "      <PointData>\n";
    os << "      </PointData>\n";
}

void VtkWriter::writeCellDataBlock(std::ostream &os)
{
    int n;
    
    os << "      <CellData>\n";
    
    n = _cellDataScalar.size();
    for (int i=0; i<n; ++i)
        writeScalarDataArrayBlock(os, _cellDataScalar[i], 8);
    
    n = _cellDataVector.size();
    for (int i=0; i<n; ++i)
        writeVectorDataArrayBlock(os, _cellDataVector[i], 8);

    os << "      </CellData>\n";
}

VtkWriterImageData::VtkWriterImageData()
{
    _x0=0; _y0=0; _z0=0;
    _dx=0; _dy=0; _dz=0;
    _nx=0; _ny=0; _nz=0;
}

void VtkWriterImageData::write(std::ostream &os) {
    int x1=0;
    int y1=0;
    int z1=0;
    int x2=x1+_nx;
    int y2=y1+_ny;
    int z2=z1+_nz;

    os << "<?xml version=\"1.0\"?>\n";
    os << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"" << _endianness << "\">\n";
    os << "  <ImageData WholeExtent=\"" << x1 << " " << x2 << " " << y1 << " " << y2 << " " << z1 << " " << z2 << "\" " \
                        "Origin=\"" << _x0 << " " << _y0 << " " << _z0 << "\" " \
                        "Spacing=\"" << _dx << " " << _dy << " " << _dz << "\">\n";
    os << "    <Piece Extent=\"" << x1 << " " << x2 << " " << y1 << " " << y2 << " " << z1 << " " << z2 << "\">\n";
    writePointDataBlock(os);
    writeCellDataBlock(os);
    os << "    </Piece>\n";
    os << "  </ImageData>\n";
    os << "</VTKFile>\n";
}

}
