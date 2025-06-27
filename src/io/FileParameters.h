/*
 * Author: ignacio.martin@imdea.org
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

#ifndef GENERALFILEPARAMETERS_H
#define GENERALFILEPARAMETERS_H

#include "Parameters.h"

namespace IO {

class FileParameters : public Parameters
{
 public:
   FileParameters(const std::string &szPath);
   FileParameters(std::istream &is) : Parameters(is) {}
   ~FileParameters();
   void dump(std::ostream &aOut, std::string const& aFilename);

  private:
   void openDir (const std::string &szPath,     const std::string &key);
   void openFile(const std::string &szFileName, const std::string &key);
};


}

#endif
