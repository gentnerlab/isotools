// $Id: StringUtils.h,v 1.4 2008/06/13 15:35:29 samn Exp $ 
#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>
#include <vector>

//remove newlines from end of str
bool Chop(char* str);

//split string based on chars in strDelim and return tokens in vstr
void Split(std::string& str,std::string& strDelim,std::vector<std::string>& vstr);
//split string based on chars in strDelim and return tokens in vstr , will modify cStr!!!!
void Split(char* cStr,std::string& strDelim,std::vector<std::string>& vstr);

//remove double quotes from str
void StripQuotes(std::string& str);

#endif
