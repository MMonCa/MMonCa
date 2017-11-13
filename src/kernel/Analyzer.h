/*
 * Analyzer.h
 *
 *  Created on: Mar 14, 2011
 *
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

#ifndef ANALYZER_H_
#define ANALYZER_H_

#include "Coordinates.h"
#include "io/Diagnostic.h"
#include "io/ParameterManager.h"

#include <cmath>
#include <map>
#include <string>
#include <cstring>
#include <cstdlib>

namespace 	Kernel {

template <typename T>
class Analyzer
{
public:
	Analyzer(const std::string &txt);

	bool        operator()(unsigned degeneracy, T x, T y, T z, unsigned size);
	Coordinates operator()(unsigned degeneracy, T p1, T p2,    unsigned size);
	T           operator()(const std::map<std::string, double> &);
	const std::string errorTxt() const;

private:
	bool isWhite(char c) { return (c==' ' || c==9); }
	bool isDelim(char c) { return (isOper(c) || c==' ' || c==9 || c=='\r' || c==0); }
	bool isOper (char c) { return
		(c=='+' || c=='-' || c=='/' || c=='*' || c=='%' || c=='^' ||
		 c=='(' || c==')' || c=='<' || c=='>' || c=='&' || c=='|' ||
		 c=='='); }
	bool isNumber(char c) { return (c >='0' && c<='9') || c=='.'; }
	void next_token();

	T start(const std::string &);

	void assign(T &res);
	void log_or(T &res);
	void log_and(T &res);
	void moreless(T &res);
	void add(T &res);
	void multiply (T &res);
	void fn(T &res);
	void power(T &res);
	void sign(T &res);
	void parenthesis(T &res);
	void primitive(T & res);

	enum TOKENS { DELIM, VARIABLE, NUMBER, FUNCTION } _token_type;
	std::map<std::string, T (*)(T)> _functions;
	T _variables[(1+'Z'-'A')*2];
	std::string _txt;
	const char *_prog, *_init;
	char * _end, *_operator;
	char _buffer[500];
	bool _lvalue;
};


template <typename T>
Analyzer<T>::Analyzer(const std::string &txt) : _txt(txt)
{
	for(int i=0; i < (1+'Z'-'A')*2; ++i)
		_variables[i] = 0;

	_variables['Z' - 'A' + 1 + 'p'-'a'] = M_PI;
	_variables['Z' - 'A' + 1 + 'e'-'a'] = M_E;

	_functions[std::string("sin")] = std::sin;
	_functions["cos"] = std::cos;
	_functions["tan"] = std::tan;
	_functions["abs"] = std::fabs;
	_functions["ln"]  = std::log;
	_functions["log"] = std::log10;
	_functions["sqrt"]= std::sqrt;
	_functions["exp"] = std::exp;
}

// x, y, z -> position
// size    -> number of defects
template <typename T>
bool Analyzer<T>::operator()(unsigned degen, T x, T y, T z, unsigned size)
{
	T res=0;
   _variables['Z' - 'A' + 1 + 'd'-'a'] = degen;
   _variables['Z' - 'A' + 1 + 's'-'a'] = size;
   _variables['Z' - 'A' + 1 + 'x'-'a'] = x;
   _variables['Z' - 'A' + 1 + 'y'-'a'] = y;
   _variables['Z' - 'A' + 1 + 'z'-'a'] = z;

   std::vector<std::string> txts;
   IO::ParameterManager::getTokens(_txt, ';', txts);
   for(std::vector<std::string>::iterator i=txts.begin(); i!=txts.end(); ++i)
	   res = start(*i);
   return res;
}

template <typename T>  //returns a value given the input variables...
T Analyzer<T>::operator()(const std::map<std::string, double> &theMap)
{
	for(std::map<std::string, double>::const_iterator it=theMap.begin(); it!=theMap.end(); ++it)
		if(it->first.size() == 1 && std::isupper(it->first[0]) )
		{
			char c = it->first[0];
			_variables[c - 'A'] = it->second;
		}
		else
			ERRORMSG("Analyzer: Variable " << it->first << " not recognized. Must be an uppercase character");

	T res=0;
	std::vector<std::string> txts;
	IO::ParameterManager::getTokens(_txt, ';', txts);
	for(std::vector<std::string>::iterator i=txts.begin(); i!=txts.end(); ++i)
		res = start(*i);

	return res;
}

// t, u -> surface parametrization
// size    -> number of defects
template <typename T>
Kernel::Coordinates Analyzer<T>::operator()(unsigned degen, T t, T u, unsigned size)
{
   _variables['Z' - 'A' + 1 + 'd'-'a'] = degen;
   _variables['Z' - 'A' + 1 + 's'-'a'] = size;
   _variables['Z' - 'A' + 1 + 't'-'a'] = t;
   _variables['Z' - 'A' + 1 + 'u'-'a'] = u;
   _variables['X'-'A'] = 0;
   _variables['Y'-'A'] = 0;
   _variables['Z'-'A'] = 0;
   std::vector<std::string> txts;
   IO::ParameterManager::getTokens(_txt, ';', txts);
   for(std::vector<std::string>::iterator i=txts.begin(); i!=txts.end(); ++i)
   	   start(*i);
   return Coordinates(_variables['X'-'A'], _variables['Y'-'A'], _variables['Z'-'A']);
}

template <typename T>
T  Analyzer<T>::start(const std::string &expr)
{
	T res = 0;

   _lvalue = true;
   _init = _prog = expr.c_str();
   if(expr.size() > 498)
	   ERRORMSG("Analizer: Expresion '" << expr << "' is too big");

   strcpy(_buffer, _prog);
   _end = _buffer;

   next_token();
   if(!_operator[0])
	  ERRORMSG("Analyzer: No expression found");

   try
   {
	   assign(res);
   }
   catch(const std::string &err)
   {
	   std::stringstream ss;
	   ss << "Analyzer: " << err << "\n" << expr << "\n";
	   for(const char *pr = expr.c_str(); pr <_prog; ++pr)
		   ss << ' ';
	   ss << "^ here";
	   ERRORMSG(ss.str());
   }
   catch(...)
   {
	   ERRORMSG("Unknown exception in Analyzer");
   }
   if(_token_type != DELIM)
	  ERRORMSG("Analyzer: Syntax error at the end of expression: " << expr);
   return res;
}

template <typename T>
void Analyzer<T>::assign(T &res)
{
	char var = _operator[0];
	log_or(res);

	char option;
	if( (option = _operator[0]) == '=')
	{
		if(!_lvalue)
			throw std::string("Analyzer: l-value exected");
		next_token();

		assign(res);
		if(var > 'Z')
			_variables['Z' - 'A' + 1 + var -'a'] = res;
		else
			_variables[var -'A'] = res;
	}
}

// |
template <typename T>
void Analyzer<T>::log_or(T &res)
{
	log_and(res);

	char option;
	while( (option = _operator[0]) == '|')
	{
		next_token();
		T temp_res;
		log_and(temp_res);
		res = (res || temp_res ? 1 : 0);
	}
}

// |
template <typename T>
void Analyzer<T>::log_and(T &res)
{
	moreless(res);

	char option;
	while( (option = _operator[0]) == '&')
	{
		next_token();
		T temp_res;
		moreless(temp_res);
		res = (res && temp_res ? 1 : 0);
	}
}

// < >
template <typename T>
void Analyzer<T>::moreless(T &res)
{
	add(res);

	char option;
	while( (option = _operator[0]) == '<' || option == '>' )
	{
		next_token();
		T temp_res;
		add(temp_res);
		res = ( ((option == '<'? temp_res - res : res - temp_res) > 0)? 1 : 0 );
	}
}

template <typename T>
void Analyzer<T>::add(T &res)
{
	multiply(res);

	char option;
	while( (option =  _operator[0]) == '+' || option == '-')
	{
		next_token();
		T temp_res;
		multiply(temp_res);
		res = (option == '+')? res+temp_res : res-temp_res;
	}
}

template <typename T>
void Analyzer<T>::multiply (T &res)
{
	fn(res);

	char option;
	while ( (option=_operator[0]) == '*' || option == '/')
	{
		next_token();
		T temp_res;
		fn(temp_res);
		res = (option == '*'? res*temp_res : res / temp_res);
	}
}

template <typename T>
void Analyzer<T>::fn(T &res)
{
	if (_token_type == FUNCTION)
	{
		std::string temp = _operator;
		next_token();
		power(res);
		typename std::map<std::string, T (*)(T)>::iterator it = _functions.find(temp);
		if(it == _functions.end())
			throw std::string("Syntax error: Function name not found "+temp);
		res = it->second(res);
	}
	else
		power(res);
}

template <typename T>
void Analyzer<T>::power(T &res)
{
	sign(res);
	while ( _operator[0] == '^')
	{
		T temp_res;
		next_token();
		sign(temp_res);
		res = std::pow(res,temp_res);
	}
}
template <typename T>
void Analyzer<T>::sign(T &res)
{
	char sign = 0;
	if ( _token_type == DELIM && (_operator[0]=='+' || _operator[0]=='-') )
	{
		sign = _operator[0];
		next_token();
	}
	parenthesis(res);
	if (sign == '-')
		res = -res;
}

// ( )
template <typename T>
void Analyzer<T>::parenthesis(T &res)
{
	if ( (_operator[0]== '(') && (_token_type == DELIM) )
	{
		next_token();
		assign(res);
		if( _operator[0] != ')')
			throw std::string("Expected ')'");
		next_token();
	}
	else
		primitive(res);
}

// variables, numbers and functions
template <typename T>
void Analyzer<T>::primitive(T & res)
{
	switch(_token_type)
	{
	case VARIABLE:
	{
		unsigned idx;
		if(*_operator <= 'Z')
			idx = *_operator - 'A';
		else
			idx = 'Z' + 1 - 'A' + *_operator - 'a';
		res = _variables[idx];
		next_token();
		break;
	}
	case NUMBER:
		res=atof(_operator);
		next_token();
		break;
	case FUNCTION:
		fn(res);
		break;
	case DELIM:
		throw std::string("Syntax error: Too many delimiters");
	default:
		throw std::string("Invalid token type");
	}
}

template <typename T>
void Analyzer<T>::next_token(void)
{
	_token_type=DELIM;
	while(isWhite(*_prog)) ++_prog;
	_operator = _end = _buffer + (_prog - _init);
	*_end = *_prog; //it might have been set to 0 before

	if(isOper(*_prog))
	{
		if(*_prog != '=')
			_lvalue = false;
		_token_type = DELIM;
		_end++; _prog++;
	}
	else if(isalpha(*_prog))
	{
		unsigned c = 0;
		while(!isDelim(*_prog)) { c++; _end++; _prog++;}
		if(c > 1)
		{
			_token_type = FUNCTION;
			_lvalue = false;
		}
		else
			_token_type = VARIABLE;
	}
	else if (isNumber(*_prog))
	{
		_lvalue = false;
		bool flag=false;
		bool dot =false;
		while (true)
		{
			switch(*_prog)
			{
				case 'e':
				case 'E':
					if (flag)
						throw std::string("Syntax error: 'e' unexpected");
					if (*(_prog+1) == '-')
						{_end++; _prog++;}
					_end++; _prog++;
					flag=true;
				break;
				case '.':
					if(dot || flag)
						throw std::string("Syntax error, '.' unexpected");
					dot = true;
				default:
					if(isDelim(*_prog))
					{
						_token_type=NUMBER;
						return;
					}
					_end++; _prog++;
				break;
			}
		}
	}
	*_end = 0;
}

}
#endif /* Analyzer_H_ */
