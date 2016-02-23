// $Id: diag.cpp,v 1.1 2011/01/31 16:06:38 samn Exp $ 
/*
This file is part of the FIRE -- Flexible Image Retrieval System

FIRE is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

FIRE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with FIRE; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <string>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <unistd.h>
#include "diag.hpp"

using namespace std;

// common function, where else to place ??
bool CheckRuntimeDebugLevel(int level, bool set) {
  static int s_Level = -1;
  if (set)  s_Level=level;
  if (s_Level == -1) {
    char * p = getenv("DEBUG_LEVEL");
    if (p) {
      string s(p);
      cout << "DEBUG_LEVEL is " << s << endl;
      istringstream is(s);
      is >> s_Level;
      cout << "DEBUG_LEVEL set to " << s_Level << endl;
    }
  }
  if (s_Level == -1) s_Level=99;
  return level<=s_Level;
}

string GetCurrentWorkingDirectory() {
  char cwd[1024];
  getcwd(cwd, 1024);
  return ::std::string(cwd);
}

void printCmdline(uint argc, char **argv) {
  cout << argv[0];
  for(uint i=1;i<argc;++i) {
    cout << " "  << argv[i];
  }
  cout << endl;
}
