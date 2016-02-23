// $Id: Log.h,v 1.2 2011/01/08 01:15:38 samn Exp $ 
/*

isoi - This program calculates the Isolation Information (IsoI) cluster quality measures
described in the reference below. These measures were designed for clusters
obtained from neural extracellular recordings, but are applicable to an
arbitrary dataset.

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

*/

#ifndef LOG_H
#define LOG_H

#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "A2D.h"

using namespace std;

bool Write2Log(const char* cstr,...);
bool Write2LogPlain(const char* cstr,...);
bool WriteVec2Log(vector<int>& v);
bool WriteVec2Log(vector<float>& v);
bool WriteVec2Log(vector<double>& v);
bool WriteMat2Log(vector<vector<float> >& m);
bool WriteMat2Log(vector<vector<double> >& m);
bool WriteMat2Log(A2D<float>& m);
bool WriteMat2Log(A2D<int>& m);
string GetDateTime(time_t t);

//stupid log file class
struct LogF
{
	FILE* m_fp;
	LogF()
		:m_fp(0)
	{}
	FILE* Open()
	{
		m_fp=fopen("__isoi__.log","a+");
		return m_fp;
	}
	void Close()
	{	if(m_fp)
		{	fclose(m_fp);
			m_fp=0;
		}
	}
	~LogF()
	{	Close();
	}
};

#endif
