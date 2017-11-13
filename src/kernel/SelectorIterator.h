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

#ifndef _SELECTOR_IT_H_
#define _SELECTOR_IT_H_
 
 namespace Kernel {

 template <typename T_EVENT, typename T_RATE>
 class SelectorIterator
 {
 public:
	SelectorIterator() {}
	SelectorIterator(const SelectorIterator<T_EVENT, T_RATE> & );
 	const SelectorIterator<T_EVENT, T_RATE> & operator++();
	T_EVENT *operator->() const { return _sel->_events[_currentIdx]; }
	T_EVENT *operator*()  const { return _sel->_events[_currentIdx]; }
	bool     operator!=(const SelectorIterator<T_EVENT, T_RATE> &it) const
	{ return (_sel != it._sel || _currentIdx != it._currentIdx); }
	bool     operator==(const SelectorIterator<T_EVENT, T_RATE> &it) const
		{ return !operator!=(it); }
	
 private:
	const Selector<T_EVENT, T_RATE> *_sel;
	int _currentIdx;
	void begin() ;
	void end()   { _currentIdx = -1; }
	
	friend class Selector<T_EVENT, T_RATE>;
};

 template <typename T_EVENT, typename T_RATE>
 SelectorIterator<T_EVENT, T_RATE>::SelectorIterator(const SelectorIterator<T_EVENT, T_RATE> &si)
 {
	 _sel = si._sel;
	 _currentIdx = si._currentIdx;
 }

template <typename T_EVENT, typename T_RATE>
const SelectorIterator<T_EVENT, T_RATE> & SelectorIterator<T_EVENT, T_RATE>::operator++()
{
	if(++_currentIdx == _sel->_howMany) 
		end();
	else while(_sel->_tags[_currentIdx] != 0 || _sel->_events[_currentIdx] == 0)
		if(++_currentIdx == _sel->_howMany)
		{ 
			end(); 
			break; 
		}
	return *this;
}

template <typename T_EVENT, typename T_RATE>
void SelectorIterator<T_EVENT, T_RATE>::begin()
{
	_currentIdx = -1;
	operator++();
}

}
#endif
