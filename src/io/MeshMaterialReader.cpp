#include "Diagnostic.h"
#include "MeshMaterialReader.h"
#include <algorithm>
#include <fstream>


MeshMaterialReader::MeshMaterialReader(std::string const& aFilename) {  // TODO +MAX_MATERIALS
    std::ifstream in(aFilename);
    if(in.good()) {
        nlohmann::json top;
        try {
            in >> top;
            _linesX = readLines(top["linesX"]);
            _linesY = readLines(top["linesY"]);
            _linesZ = readLines(top["linesZ"]);
            _originalMap = top["materialMapping"];
            auto const materials = top["materialIDs"];
            if(materials.size() != (_linesX.size() - 1u) * (_linesY.size() - 1u) * (_linesZ.size() - 1u)) {
                throw std::invalid_argument("materialIDs.size() != (linesX.size() - 1) * (linesY.size() - 1) * (linesZ.size() - 1)");
            }
            _materials.reserve(materials.size());
            for(auto const& value : materials) {
                float raw = value;
                _materials.push_back(raw + Kernel::MAX_MATERIALS);
            }
        }
        catch(std::exception &e) {
            ERRORMSG("Error parsing " << aFilename << ": " << e.what());
        }
    }
    else {
		ERRORMSG("Could not find the mesh description file " << aFilename);
    }
}

void MeshMaterialReader::translate(std::string const& aFrom, Kernel::M_TYPE const aTo) {
    auto const from = static_cast<float>(_originalMap[aFrom]) + Kernel::MAX_MATERIALS;
    std::replace(_materials.begin(), _materials.end(), static_cast<Kernel::M_TYPE>(from), aTo);
}

std::vector<float> MeshMaterialReader::readLines(nlohmann::json const& aJson) {
    std::vector<float> result;
    for(auto const& value : aJson) {
        float raw = value;
        if(result.size() > 0u && result.back() >= raw) {
            throw std::invalid_argument("lines* array must be strictly increasing.");
        }
        result.push_back(raw);
    }
    if(result.size() < 2u) {
        throw std::invalid_argument("lines* array must have at least 2 values.");
    }
    return result;
}