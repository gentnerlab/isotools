/* $Id: MyObj.h,v 1.2 2011/01/07 04:22:01 samn Exp $ */
// MyObj.h: interface for the CMyObject class.
//
/*

isorat - This program calculates the Isolation Distance (IsoI) and L-Ratio cluster
quality measures described in the references below. These measures were
designed for clusters obtained from neural extracellular recordings.

Copyright (C) 2003-2011 Sam Neymotin & BioSignal Group

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

For more information contact Sam Neymotin ( samn at neurosim dot downstate dot edu )
 or Andre Fenton ( afenton at nyu dot edu ).

References:
 The methods used in this program are described in an article in press
  at The Journal of Neuroscience,
  Measuring the quality of neuronal identification in ensemble recordings
  by Neymotin SA, Lytton WW, Olypher AO, Fenton AA (2011).

  and in, Quantitative measures of cluster quality for use in extracellular recordings
  Schmitzer-Torbert N, Jackson J, Henze D, Harris K, Redish AD
  Neuroscience 131(1):1-11, 2005.

*/

//////////////////////////////////////////////////////////////////////

#ifndef MYOBJ_H_
#define MYOBJ_H_

#include <vector>
#include <deque>
#include <set>
#include <memory>
#include <float.h>

using namespace std;

////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
// CMyObject
class CMyObject 
{
public:
	CMyObject();
	virtual ~CMyObject();
	virtual void GetValue(){};
};

////////////////////////////////////////////////////////////////////////
// Containers
typedef vector <float> VERTEX;
typedef vector <int> MY_INT_VECT;
typedef deque <CMyObject*> MY_STACK;
typedef deque <int> MY_INT_STACK;

//probability type
typedef float prob_t;

//#if sizeof(prob_t) = sizeof(float)
  const float MAX_PROB_T = FLT_MAX;
  const float MIN_PROB_T = FLT_MIN;
//#else
//  const double MAX_PROB_T = DBL_MAX;
//  const double MIN_PROB_T = DBL_MIN;
//#endif

const prob_t one = 1.0, two = 2.0;

#endif 

