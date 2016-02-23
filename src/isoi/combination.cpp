// $Id: combination.cpp,v 1.2 2011/01/07 04:21:12 samn Exp $ 

#include "combination.h"
#include "Log.h"

#include <iostream>
using namespace std;

using namespace btb;

int comb_main()
{
	LogF F; FILE* fp = F.Open();
  int a[7] = {1, 2, 3, 4, 5, 6, 7};

  int c = 0;
  int len = 5;
  do
    {
      fprintf(fp,"--[\t");
      for (int i = 0; i < c; ++i)
	fprintf(fp,"%d\t",a[i]);
	  fprintf(fp,"   <>\t");
      for (int i = c; i < len; ++i)
		fprintf(fp,"%d\t",a[i]);
      fprintf(fp,"]\n");
    }
  while (next_combination(a, a+c, a+len) || (++c <= len));

  c = 0;
  len = 5;
  a[1] = 1;
  a[4] = 4;
  do
    {
      fprintf(fp,"--[\t");
      for (int i = 0; i < c; ++i)
	fprintf(fp,"%d\t",a[i]);
      fprintf(fp,"   <>\t");
      for (int i = c; i < len; ++i)
	fprintf(fp,"%d\t",a[i]);
      fprintf(fp,"]\n");
    }
  while (next_combination(a, a+c, a+len) || (++c <= len));
  return 0;
}

struct jnkjnk
{
	jnkjnk()
	{comb_main();
	}
} ;//gggggg;
