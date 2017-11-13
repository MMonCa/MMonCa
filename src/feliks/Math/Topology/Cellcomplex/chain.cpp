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
/*
 *  chain.cpp
 *  feliks
 *
 *  Created by Ignacio Romero on 10/9/11.
 *  Copyright 2011 Universidad Politécnica de Madrid. All rights reserved.
 *
 */

#include "chain.h"
#include "cell.h"

using namespace feliks::topology;

baseChain::baseChain() :
	dimension(-1)
{
}

baseChain::baseChain(const int d) :
dimension(d)
{
}


baseChain::~baseChain()
{
	(*this).clear();
}


baseChain::baseChain(const baseChain& bs) :
	dimension(bs.dimension)
{
	std::map<cell *,int>::const_iterator iter=bs.begin();
	while (iter != bs.end())
	{
		(*this)[iter->first] = iter->second;
		++iter;
	}
	
}



baseChain& baseChain::operator=(const baseChain& rhs)
{
	this->clear();
	this->dimension = rhs.dimension;
	
	std::map<cell *,int>::const_iterator iter=rhs.begin();
	while (iter != rhs.end())
	{
		(*this)[iter->first] = iter->second;
		++iter;
	}
	
	return *this;
}



bool baseChain :: isZero() const
{
	int sum=0;
	
	std::map<cell *,int>::const_iterator iter=begin();
	while (iter != end())
	{
		sum += iter->second * iter->second;
		++iter;
	}
	return sum == 0;
}



void baseChain::print(std::ostream& of) const
{
	of << "\n Cells:";
	
	std::map<cell *,int>::const_iterator iter = begin();
	while ( iter != end() )
	{
        of << "\n\t";
		(*iter).first->print(of);
		of << ", with multiplicity " << (*iter).second;
		++iter;
	}
}



chain::chain() :
	baseChain()
{
}



chain::chain(const int d) :
	baseChain(d)
{
}



int chain::operator()(const cochain& co)
{
	int sum = 0;
	std::map<cell *,int>::const_iterator iter=co.begin();
	while ( iter != co.end() )
	{
		sum += (*this)[iter->first] * iter->second;
		++iter;
	}
	return sum;
}


chain chain::boundary()
{
	chain ch(dimension-1);
	
	std::map<cell *,int>::const_iterator iter = (*this).begin();
	while( iter != (*this).end() )
	{
		ch += (*iter).first->boundary() * (*iter).second;
		++iter;
	}
	
	return ch;
}


chain chain::operator+(const chain &right)
{
	chain  ch =(*this);
	return ch += right;
}




chain& chain::operator+=(const chain &right)
{	
	std::map<cell *,int>::const_iterator iter=right.begin();
	while( iter != right.end() )
	{
		(*this)[ iter->first ] += iter->second;
		++iter;
	}
	
	return *this;
}





chain chain::operator-(const chain &right)
{
	chain  ch = (*this);
	return ch -= right;
}


chain& chain::operator-=(const chain &right)
{	
	std::map<cell *,int>::const_iterator iter=right.begin();
	while( iter != right.end() )
	{
		(*this)[ iter->first ] -= iter->second;
		++iter;
	}
	
	return *this;
}


chain& chain::operator*=(const int n)
{	
	std::map<cell *,int>::const_iterator iter=this->begin();
	while( iter != this->end() )
	{
		(*this)[ iter->first ] *= n;
		++iter;
	}
	
	return *this;
}



chain chain::operator*(const int a)
{
	chain ch = (*this);
	return ch *= a;
}



void chain::print(std::ostream& of) const
{
	of << "\n Chain of size " << size();
	baseChain::print(of);
}	



cochain::cochain() :
	baseChain()
{
}


cochain::cochain(const int d) :
	baseChain(d)
{
}





cochain cochain::coboundary()
{
	cochain cobo(dimension+1);
	
	std::map<cell *,int>::iterator iter = (*this).begin();
	while( iter != (*this).end() )
	{
		cobo += (*iter).first->coboundary() * (*iter).second;
		++iter;
	}
	
	return cobo;
}




int cochain::operator()(const chain& co)
{
	int sum = 0;
	std::map<cell *,int>::const_iterator iter=co.begin();
	while ( iter != co.end() )
	{
		sum += (*this)[iter->first] * iter->second;
		++iter;
	}
	return sum;
}




cochain& cochain::operator+=(const cochain &right)
{	
	std::map<cell *,int>::const_iterator iter=right.begin();
	while( iter != right.end() )
	{
		(*this)[ iter->first ] += iter->second;
		++iter;
	}
	
	return *this;
}



cochain cochain::operator+(const cochain &right)
{
	cochain ch =(*this);
	return  ch += right;
}





cochain cochain::operator-(const cochain &right)
{
	cochain ch = (*this);
	return  ch -= right;
}


cochain& cochain::operator-=(const cochain &right)
{	
	std::map<cell *,int>::const_iterator iter=right.begin();
	while( iter != right.end() )
	{
		(*this)[ iter->first ] -= iter->second;
		++iter;
	}
	
	return *this;
}


cochain& cochain::operator*=(const int n)
{	
	std::map<cell *,int>::const_iterator iter=this->begin();
	while( iter != this->end() )
	{
		(*this)[ iter->first ] *= n;
		++iter;
	}
	
	return *this;
}



cochain cochain::operator*(const int a)
{
	cochain ch;
	ch = (*this);
	ch *= a;
	return ch;	
}



void cochain::print(std::ostream& of) const
{
	of << "\n Cochain of size " << size();
	baseChain::print(of);
}	



