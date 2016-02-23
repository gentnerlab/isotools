// $Id: StringUtils.cpp,v 1.3 2008/06/13 15:35:16 samn Exp $ 
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

#include "StringUtils.h"
#include <string.h>

using namespace std;

//remove newlines from end of str
bool Chop(char* str)
{
	int iLen = strlen(str);
	if(iLen)
	{
		if(str[iLen-1]=='\n' || str[iLen-1]=='\r')
		{
			str[iLen-1]=0;
			Chop(str); //in case \r\n
			return true;
		}
	}
	return false;
}

//split string based on chars in strDelim and return tokens in vstr
void Split(string& str,string& strDelim,vector<string>& vstr)
{
	vstr.resize(0);
	char* cStr = new char[str.length()+1];
	strcpy(cStr,str.c_str());

	char* token = strtok( cStr, strDelim.c_str() );
	while( token != NULL )
	{
		vstr.push_back(string(token));
		token = strtok( NULL, strDelim.c_str() );
	}
	delete [] cStr;
}

//split string based on chars in strDelim and return tokens in vstr
void Split(char* cStr,string& strDelim,vector<string>& vstr)
{
	vstr.resize(0);
	char* token = strtok( cStr, strDelim.c_str() );
	while( token != NULL )
	{
		vstr.push_back(string(token));
		token = strtok( NULL, strDelim.c_str() );
	}
}

void StripQuotes(string& str)
{
	string strTmp;
	int iSz = str.size(), i = 0;
	for(;i<iSz;i++)
		if(str[i]!='\"')
			strTmp+=str[i];
	str=strTmp;
}

