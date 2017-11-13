/* Author: ignacio.martin@imdea.org
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

#ifndef _SELECTOR__H_
#define _SELECTOR__H_

#include <vector>
#include <cmath>
#include <cassert>
#include "io/Diagnostic.h"
#include <set>

namespace Kernel
{
  template <typename T_EVENT, typename T_RATE>
  class SelectorIterator;

  template <typename T_EVENT, typename T_RATE>
  class Selector
  {
  public:
    typedef SelectorIterator<T_EVENT,T_RATE> iterator;
    
    Selector(int starting_levels = 1000);
    ~Selector() {};  
    
    T_RATE getTotalRate() const { return _values[0][0]; }
    void   getValues(unsigned idx, T_EVENT *&, int &tag, T_RATE &val) const;
    int insert(T_EVENT *, T_RATE, int tag);
    void remove(int idx);
    void update(int idx, T_RATE);
    void update(); //all
    void select(T_RATE rand, T_EVENT *&, int &tag) const;
    
    iterator begin() const;
    iterator end()   const;

  private:
    void increase();

    std::vector<std::vector<T_RATE> > _values;
    std::vector<T_EVENT *> _events;
    std::vector<int>       _tags;
    std::vector<int>       _holes;
    int _levels;
    int _howMany;
    
    friend class SelectorIterator<T_EVENT,T_RATE>;
  };

  template <typename T_EVENT, typename T_RATE>
  Selector<T_EVENT, T_RATE>::Selector(int howmany)
  {
    _levels = int(log2(howmany)+2);
    _values.resize(_levels);
 
    _howMany=1;
    for(int i=0; i<_levels; ++i)
    {
      _values[i].resize(_howMany, 0);
      _howMany *= 2;
    }
    _howMany /= 2;
    _events.resize(_howMany, 0);
    _holes.resize(_howMany, 0);
    _tags.resize(_howMany, 0);
    for(int i=0; i<_howMany; ++i)
      _holes[i]=i;
  }

  template <typename T_EVENT, typename T_RATE>
  int Selector<T_EVENT, T_RATE>::insert(T_EVENT *eve, T_RATE rate, int tag)
  {
    if(!_holes.size())
      increase();
    int idx = _holes.back();
    _holes.pop_back();
    _tags[idx] = tag;
    _events[idx] = eve;
    _values[_levels-1][idx] = rate;
    //update tree
    int half=idx;
    for(int lev=_levels-2; lev >=0; --lev)
    {
      half >>= 1; // /= 2;
      _values[lev][half] = _values[lev+1][half << 1] + _values[lev+1][(half << 1)+1];
    }
    return idx;
  }

  template <typename T_EVENT, typename T_RATE>
  void Selector<T_EVENT, T_RATE>::remove(int idx)
  {
    assert(_events[idx] != 0); //Removing and event that was not inserted
    _holes.push_back(idx);
    _events[idx] = 0;
    _values[_levels-1][idx] = 0;
    //update tree
    int half=idx;
    for(int lev=_levels-2; lev >=0; --lev)
    {
      half /= 2;
      _values[lev][half] = _values[lev+1][half*2] + _values[lev+1][half*2+1];
    }
  }

  template <typename T_EVENT, typename T_RATE>
  void Selector<T_EVENT, T_RATE>::select(T_RATE rand, T_EVENT *&event, int &tag) const
  {
    T_RATE pickup = rand*_values[0][0];
    int idx=0;
    int lev=1;
    for(; lev <_levels-1; ++lev)
    {
      T_RATE prev = _values[lev][idx];
      if((pickup < prev && prev != 0) || (_values[lev][idx+1] == 0)) //last one due to rounding errors
        idx *=2;
      else
      {
        idx++;
        idx *=2;
        pickup -= prev;
      }
    }
    if (pickup >=  _values[lev][idx] && _values[lev][idx+1] != 0)
      idx++;
    tag = _tags[idx];
    event = _events[idx];
  }

  template <typename T_EVENT, typename T_RATE>
  void Selector<T_EVENT, T_RATE>::increase()
  {
    _values.push_back(std::vector<T_RATE>());
    _howMany *= 2;
    //MEDMSG("Relocating space for " << _howMany << " elements");
    _values[_levels].resize(_howMany,0);
    _tags.resize(_howMany,-1);
    _events.resize(_howMany,0);
    _levels++;
    int i=0;
    for(; i<_howMany/2; ++i)
      _values[_levels-1][i] = _values[_levels-2][i];
    for(; i<_howMany; ++i)
      _holes.push_back(i);
    update();
  }

  template <typename T_EVENT, typename T_RATE>
  void Selector<T_EVENT, T_RATE>::update()
  {
    for(int i=_levels-2; i>=0; i--)
      for(unsigned j=0; j<_values[i].size(); ++j)
        _values[i][j] = _values[i+1][j*2] + _values[i+1][j*2+1];
  }
  
  template <typename T_EVENT, typename T_RATE>
  void Selector<T_EVENT, T_RATE>::update(int idx, T_RATE rate)
  {
    assert(_values[_levels-1].size() > idx);
    _values[_levels-1][idx] = rate;
    //update tree
    int half=idx;
    for(int lev=_levels-2; lev >=0; --lev)
    {
      half /= 2;
      _values[lev][half] = _values[lev+1][half*2] + _values[lev+1][half*2+1];
    }
  }
  
  template <typename T_EVENT, typename T_RATE>
  void Selector<T_EVENT, T_RATE>::getValues(unsigned idx, T_EVENT *&ev, int &tag, T_RATE &val) const
  {
  	ev = _events[idx];
	tag = _tags[idx];
	val = _values[_levels-1][idx];
  }

  template <typename T_EVENT, typename T_RATE>
  SelectorIterator<T_EVENT, T_RATE> Selector<T_EVENT, T_RATE>::begin() const
  {
  	iterator it;
	it._sel = this;
	it.begin();
	return it;
  }
  
  template <typename T_EVENT, typename T_RATE>
  SelectorIterator<T_EVENT, T_RATE> Selector<T_EVENT, T_RATE>::end() const
  {
  	iterator it;
	it._sel = this;
	it.end();
	return it;
  }
}
#include "SelectorIterator.h"
#endif
