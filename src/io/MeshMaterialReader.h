#ifndef MESHMATERIALREADER_H_
#define MESHMATERIALREADER_H_

#include "kernel/Material.h"
#include "thirdparty/json.hpp"
#include <cstdint>
#include <string>
#include <vector>


class MeshMaterialReader final {
public:
    MeshMaterialReader(std::string const& aFilename);

    std::vector<float> const& getLinesX() const { return _linesX; }
    std::vector<float> const& getLinesY() const { return _linesY; }
    std::vector<float> const& getLinesZ() const { return _linesZ; }

    void translate(std::string const& aFrom, Kernel::M_TYPE const aTo);

    Kernel::M_TYPE getMaterial(uint32_t const aIndex) const { return _materials[aIndex]; }
    std::vector<Kernel::M_TYPE> const& getMaterials() const { return _materials; }

private:
    static std::vector<float> readLines(nlohmann::json const& aJson);

    std::vector<float>          _linesX;
    std::vector<float>          _linesY;
    std::vector<float>          _linesZ;
    std::vector<Kernel::M_TYPE> _materials;
    nlohmann::json              _originalMap;
};

#endif

/*
>>   {
>>      "linesX": [
>>          0.0,
>>          1.0,
>>          2.0
>>      ],
>>      "linesY": [
>>          2.0,
>>          4.0,
>>          6.0
>>      ],
>>      "linesZ": [
>>          2.0,
>>          3.0,
>>          4.0
>>      ],
>>      "materialIDs": [
>>          1,
>>          1,
>>          1,
>>          1,
>>          0,
>>          0,
>>          2,
>>          2
>>      ],
>>      "materialMapping": {
>>          "Silicon": 0,
>>          "SiO2": 1,
>>          "SiliconGermanium": 2
>>      }
>>   }
*/