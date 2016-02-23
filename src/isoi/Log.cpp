// $Id: Log.cpp,v 1.3 2011/01/08 01:15:32 samn Exp $ 
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

#include <time.h>
#include <string>
#include <stdio.h>
#include <stdarg.h>
#include "Log.h"

using namespace std;

string GetDateTime(time_t t)
{
  struct tm	*date;		/* Date/time value */
  static time_t	last_time = -1;	/* Last time value */
  static char	s[1024];	/* Date/time string */
  static const char * const months[12] =
		{		/* Months */
		  "Jan",
		  "Feb",
		  "Mar",
		  "Apr",
		  "May",
		  "Jun",
		  "Jul",
		  "Aug",
		  "Sep",
		  "Oct",
		  "Nov",
		  "Dec"
		};


  if (t != last_time)
  {
    last_time = t;
  
    date = localtime(&t);

    sprintf(s,  "%02d/%s/%04d:%02d:%02d:%02d",
	     date->tm_mday, months[date->tm_mon], 1900 + date->tm_year,
	     date->tm_hour, date->tm_min, date->tm_sec);
  }

  return string(s);
}

bool Write2LogPlain(const char* cstr,...)
{
	try
	{		
		const char* fname = "__isoi__.log";
		FILE* fp = fopen(fname,"a+");
		if(!fp)
		  {	fprintf(stderr,"Couldn't write to __isoi__.log!!!\n");
			return false;
		}

		va_list ap;

		static char line[32768]={0};

		va_start(ap, cstr);
		vsprintf(line, cstr, ap);
		va_end(ap);

		fprintf(fp,"%s",line);

		fclose(fp);

		return true;
	}
	catch(...)
	  {	fprintf(stderr,"Exception in Write2Log!!!\n");
		return false;
	}
}

bool Write2Log(const char* cstr,...)
{
	try
	{		
		const char* fname = "__isoi__.log";
		FILE* fp = fopen(fname,"a+");
		if(!fp)
		  {	fprintf(stderr,"Couldn't write to __isoi__.log!!!\n");
			return false;
		}

		va_list ap;

		static char line[32768]={0};

		va_start(ap, cstr);
		vsprintf(line, cstr, ap);
		va_end(ap);

		fprintf(fp,"[%s]: %s\n",GetDateTime(time(0)).c_str(),line);

		fclose(fp);

		return true;
	}
	catch(...)
	  {	fprintf(stderr,"Exception in Write2Log!!!\n");
		return false;
	}
}


bool WriteVec2Log(vector<float>& v)
{
	try
	{
		string str("\n"); char tmp[1024];
		int iSz = v.size() , i;
		sprintf(tmp,"vsz=%d\n",iSz);
		str+=tmp;
		for(i=0;i<iSz;i++)
		{
			sprintf(tmp,"%.4f\t",v[i]);
			str += tmp;
		}
		return Write2Log(str.c_str());
	}
	catch(...)
	{
		return false;
	}
}

bool WriteVec2Log(vector<double>& v)
{
	try
	{
		std::string str("\n");
		char tmp[1024];
		int iSz = v.size() , i;
		sprintf(tmp,"vsz=%d\n",iSz);
		str+=tmp;
		Write2Log(str.c_str());

		const char* fname = "__isoi__.log";
		FILE* fp = fopen(fname,"a+");
		if(!fp)
			return false;

		for(i=0;i<iSz;i++)
			fprintf(fp,"%.4f\t",v[i]);
		fprintf(fp,"\n");

		fclose(fp);
		
		return true;
	}
	catch(...)
	{	try{Write2Log("Exception in WriteVec2Log!");}catch(...){}
		return false;
	}
}

bool WriteMat2Log(vector<vector<float> >& m)
{
	try
	{
		string str("\n"); char tmp[1024];
		int iRows = m.size(), iCols = m[0].size() , i, j;
		sprintf(tmp,"%d X %d\n",iRows,iCols);
		str+=tmp;
		for(i=0;i<iRows;i++)
		{	
			for(j=0;j<iCols;j++)
			{
				sprintf(tmp,"%.4f\t",m[i][j]);
				str += tmp;
			}
			str += "\n";
		}
		return Write2Log(str.c_str());
	}
	catch(...)
	{
		return false;
	}
}

bool WriteMat2Log(A2D<float>& m)
{
	try
	{
		string str("\n"); char tmp[1024];
		int iRows = m.Rows(), iCols = m.Cols() , i, j;
		sprintf(tmp,"%d X %d\n",iRows,iCols);
		str+=tmp;
		for(i=0;i<iRows;i++)
		{	
			for(j=0;j<iCols;j++)
			{
				sprintf(tmp,"%.4f\t",m[i][j]);
				str += tmp;
			}
			str += "\n";
		}
		return Write2Log(str.c_str());
	}
	catch(...)
	{
		return false;
	}
}

bool WriteMat2Log(A2D<int>& m)
{
	try
	{
		string str("\n"); char tmp[1024];
		int iRows = m.Rows(), iCols = m.Cols() , i, j;
		sprintf(tmp,"%d X %d\n",iRows,iCols);
		str+=tmp;
		for(i=0;i<iRows;i++)
		{	
			for(j=0;j<iCols;j++)
			{
				sprintf(tmp,"%d\t",m[i][j]);
				str += tmp;
			}
			str += "\n";
		}
		return Write2Log(str.c_str());
	}
	catch(...)
	{
		return false;
	}
}


bool WriteMat2Log(vector<vector<double> >& m)
{
	try
	{
		string str("\n"); char tmp[1024];
		int iRows = m.size(), iCols = m[0].size() , i, j;
		sprintf(tmp,"%d X %d\n",iRows,iCols);
		str+=tmp;
		for(i=0;i<iRows;i++)
		{	
			for(j=0;j<iCols;j++)
			{
				sprintf(tmp,"%.4f\t",m[i][j]);
				str += tmp;
			}
			str += "\n";
		}
		return Write2Log(str.c_str());
	}
	catch(...)
	{
		return false;
	}
}

bool WriteVec2Log(vector<int>& v)
{
	try
	{
		std::string str("\n");
		char tmp[1024];
		int iSz = v.size() , i;
		sprintf(tmp,"vsz=%d\n",iSz);
		str+=tmp;
		Write2Log(str.c_str());

		const char* fname = "__isoi__.log";
		FILE* fp = fopen(fname,"a+");
		if(!fp)
			return false;

		for(i=0;i<iSz;i++)
			fprintf(fp,"%d\t",v[i]);
		fprintf(fp,"\n");

		fclose(fp);
		
		return true;
	}
	catch(...)
	{	try{Write2Log("Exception in WriteVec2Log!");}catch(...){}
		return false;
	}
}
