// $Id: StringUtils.cpp,v 1.3 2008/06/13 15:35:16 samn Exp $ 

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

