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

#ifndef GENERALPARAMETERS_H
#define GENERALPARAMETERS_H
#include <string>
#include <vector>
#include <map>
#include "Diagnostic.h"
#include "ArrheniusAlloys.h"
#include "Polynomial.h"
#include "kernel/Coordinates.h"

struct Tcl_Interp;

namespace IO
{    
	template <typename S, typename T>
	class array : public std::vector<std::pair<S, T> >
	{
	};

class Parameters
{
public:
	enum AATYPE { AA_FULL, AA_PREFACTOR, AA_ENERGY };
	Parameters() {}
	Parameters(std::istream &);
	virtual ~Parameters();

    float                                      getFloat(const std::string &key) const;
    std::vector<float>                         getFloats(const std::string &key) const;
    int                                  	   getInt(const std::string &key) const;
    std::string                          	   getString(const std::string &key) const;
    std::vector<std::string>				   getStrings(const std::string &key, unsigned howMany) const;
    std::map<std::string, std::string>   	   getStringMap(const std::string &key) const;
    std::map<std::string, IO::ArrheniusAlloys> getArrheniusAlloysMap(const std::string &key) const;
    std::map<std::string, float>               getFloatMap(const std::string &key) const;
    std::map<std::string, int>                 getIntMap(const std::string &key) const;
    std::map<std::string, bool>                getBoolMap(const std::string &key) const;
    array<std::string,std::string>             getArray(const std::string &key) const;
    bool                                       getBool(const std::string &key) const;
    bool 			                           specified(const std::string &key, const std::string &type) const;
    bool 								       specified(const std::string &key) const;
    ArrheniusAlloys                            getArrheniusAlloys(const std::string &key) const;
    Kernel::Coordinates                        getCoordinates(const std::string &key) const;
    Polynomial                       		   getPolynomial(const std::string &key) const;

    void                                       setCommand  (const std::string &type, const std::string &key, const std::string &value, const std::string &index);
    void                                       addCommand  (const std::string &type, const std::string &key, const std::string &value, const std::string &index);
    std::string              				   getCommand  (const std::string &type, const std::string &key, const std::string &index) const;
    void								       unsetCommand(const std::string &type, const std::string &key, const std::string &index);

    /* load Procedure requires the number of parameters in the procedure. The rest require an "idx". For several parameters, idx is
     * a string with the parameters separated with spaces: arg0 arg1 arg2 ... No check is done here, the check is done by Tcl.
     */
    void			                           loadProcedure   		  (Tcl_Interp *pTcl, const std::string &key, unsigned argc) const;

    std::map<std::string, ArrheniusAlloys>     getArrheniusAlloysProc (Tcl_Interp *pTcl, const std::string &key, AATYPE type) const;
    std::map<std::string, Arrhenius>           getArrheniusProc       (Tcl_Interp *pTcl, const std::string &key) const;
    float        						       getFloatProc    		  (Tcl_Interp *pTcl, const std::string &key, const std::string &idx) const;
    float									   getFloatProcArgs       (Tcl_Interp *pTcl, const std::string &key, const std::string &args) const;
    std::map<std::string, float>			   getFloatProc    		  (Tcl_Interp *pTcl, const std::string &key) const;

    static ArrheniusAlloys				       toArrheniusAlloys(const std::string &key, const std::string &txt, AATYPE);
    static std::string                         getTextInBrackets(std::stringstream &);

    void restart(std::ostream &) const;
protected:
    // key -> (type, value)
    typedef std::map<std::string, std::pair<std::string, std::string> > mmap;
    mmap _map;
    mutable std::map<std::string, bool> _used; //to keep track if they are used.
    void insert(const std::string &key, const std::string &type, const std::string &value);

private:
    const std::string &                get(const std::string &key, const std::string &type) const;

};

}

#endif
