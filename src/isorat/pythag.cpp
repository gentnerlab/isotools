// $Id: pythag.cpp,v 1.1 2008/05/06 02:29:30 samn Exp $ 

#include <cmath>
#include "nr.h"
using namespace std;

DP NR::pythag(const DP a, const DP b)
{
	DP absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

DP NR::pythag(const float a, const float b)
{
	float absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
