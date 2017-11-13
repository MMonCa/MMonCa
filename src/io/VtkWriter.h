/**
* \file VtkWriter.h
* \brief Vtk writer interface.
* \date December 01, 2014
* \author Benoit SKlenard
*/

#ifndef VTK_WRITER_H
#define VTK_WRITER_H

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <stdint.h>

namespace IO {

class VtkWriter {
    public:
        enum VtkFormat { ASCII, BINARY };
        enum VtkType { INT8, UINT8, INT16, UINT16, INT32, UINT32, INT64, UINT64, FLOAT32, FLOAT64 };

        VtkWriter();
        virtual ~VtkWriter() {}

        void setFileName(const std::string &filename) { _filename = filename; }
        
        void write(); // = 0;
        virtual void write(std::ostream &os) {} // = 0;

        inline void addCellData(const std::string &name, const std::vector<int8_t> &data, VtkFormat format=BINARY);
        inline void addCellData(const std::string &name, const std::vector<uint8_t> &data, VtkFormat format=BINARY);
        inline void addCellData(const std::string &name, const std::vector<int16_t> &data, VtkFormat format=BINARY);
        inline void addCellData(const std::string &name, const std::vector<uint16_t> &data, VtkFormat format=BINARY);
        inline void addCellData(const std::string &name, const std::vector<int32_t> &data, VtkFormat format=BINARY);
        inline void addCellData(const std::string &name, const std::vector<uint32_t> &data, VtkFormat format=BINARY);
        inline void addCellData(const std::string &name, const std::vector<int64_t> &data, VtkFormat format=BINARY);
        inline void addCellData(const std::string &name, const std::vector<uint64_t> &data, VtkFormat format=BINARY);
        inline void addCellData(const std::string &name, const std::vector<float> &data, VtkFormat format=BINARY);
        inline void addCellData(const std::string &name, const std::vector<double> &data, VtkFormat format=BINARY);

        inline void addCellData(const std::string &name, const std::vector<std::vector<float> > &data, int nComponents, VtkFormat format=ASCII);
        inline void addCellData(const std::string &name, const std::vector<std::vector<double> > &data, int nComponents, VtkFormat format=ASCII);

    protected:
        struct DataSet {
            std::string _name;
            VtkFormat _format;
            VtkType _type;
            void *_data;
            int _nComponents;
        };

        std::vector<DataSet> _cellDataScalar;
        std::vector<DataSet> _cellDataVector;
        std::string _filename;
        std::string _endianness;

        static const std::string vtkTypeString[];
        static const std::string vtkFormatString[];
       
        void writePointDataBlock(std::ostream &os);
        void writeCellDataBlock(std::ostream &os);
        
        void writeScalarDataArrayBlock(std::ostream &os, const DataSet &set, int shift=0);
        void writeVectorDataArrayBlock(std::ostream &os, const DataSet &set, int shift=0);
        
        template <typename T> void writeBinaryScalars(std::ostream &os, const std::vector<T> &v);
        
        template <typename T> void writeAsciiScalars(std::ostream &os, const std::vector<T> &v, int shift);
        template <typename T> void writeAsciiVectors(std::ostream &os, const std::vector<std::vector<T> > &v, int nComponents, int shift);
        
        void indent(std::ostream &os, int shift);
};

inline void VtkWriter::addCellData(const std::string &name, const std::vector<int8_t> &data, VtkFormat format) {
    DataSet d = { name, format, INT8, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<uint8_t> &data, VtkFormat format) {
    DataSet d = { name, format, UINT8, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<int16_t> &data, VtkFormat format) {
    DataSet d = { name, format, INT16, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<uint16_t> &data, VtkFormat format) {
    DataSet d = { name, format, UINT16, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<int32_t> &data, VtkFormat format) {
    DataSet d = { name, format, INT32, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<uint32_t> &data, VtkFormat format) {
    DataSet d = { name, format, UINT32, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<int64_t> &data, VtkFormat format) {
    DataSet d = { name, format, INT64, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<uint64_t> &data, VtkFormat format) {
    DataSet d = { name, format, UINT64, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<float> &data, VtkFormat format) {
    DataSet d = { name, format, FLOAT32, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<double> &data, VtkFormat format) {
    DataSet d = { name, format, FLOAT64, (void *)&data, 0};
    _cellDataScalar.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<std::vector<float> > &data, int nComponents, VtkFormat format) {
    DataSet d = { name, format, FLOAT32, (void *)&data, nComponents};
    _cellDataVector.push_back(d);
}

inline void VtkWriter::addCellData(const std::string &name, const std::vector<std::vector<double> > &data, int nComponents, VtkFormat format) {
    DataSet d = { name, format, FLOAT64, (void *)&data, nComponents};
    _cellDataVector.push_back(d);
}


class VtkWriterImageData : public VtkWriter {
    public:
        VtkWriterImageData();
        void setStructure(double x0, double y0, double z0, double dx, double dy, double dz, int nx, int ny, int nz);
        
        void setOrigin(double x0, double y0, double z0) { _x0=x0; _y0=y0; _z0=z0; }
        void setSize(int nx, int ny, int nz) { _nx=nx; _ny=ny; _nz=nz; }
        void setSpacing(double dx, double dy, double dz) { _dx=dx; _dy=dy; _dz=dz; }

        void write();
        void write(std::ostream &os);

    private:
        int _nx, _ny, _nz; // number of lines
        
        double _x0, _y0, _z0; // origin
        double _dx, _dy, _dz; // spacing
};

class VtkWriterRectilinearGrid : public VtkWriter {
    public:
        VtkWriterRectilinearGrid() {} 
    private:
        std::vector<double> _xlines;
        std::vector<double> _ylines;
        std::vector<double> _zlines;
};

class VtkWriterStructuredGrid : public VtkWriter { 
    public:
        VtkWriterStructuredGrid() {}
    
    private:
};

class VtkWriterUnstructuredGrid : public VtkWriter {
    public:
        VtkWriterUnstructuredGrid() {}

    private:
};

}

#endif
