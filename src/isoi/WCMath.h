// $Id: WCMath.h,v 1.2 2011/01/08 01:13:24 samn Exp $ 
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

*/

#pragma once
#include <cmath>
#include <vector>
#include "A2D.h"

using namespace std;

inline long Factorial(long n)
{
	if(n <= 1) return 1;
	return n * Factorial(n - 1);
}

inline long DFactorial(long n)
{
	if(n<=1) return 1;
	return n * DFactorial(n-2);
}

//NR version of logarithm of gamma function
inline float lngamma(float xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

inline double Gamma(double R)
{
	int iR = (int) R;
	double d = iR;
	if(d == R)
	{
		return Factorial(iR-1);
	}
	else
	{
		const double PI=3.14159265358979323846;
		//R = n/2 + 1;
		//n = (R - 1)*2;		
		int n = static_cast<int>(( R - 1.0 ) * 2.0);
		return (sqrt(PI) * DFactorial(n)) / pow(2,(n+1)/2);
	}
}

//monotonically decreasing
template< class T >
inline bool MonoDec(std::vector<T>& v)
{
	int isz = v.size(), i, j;
	for(i=0;i<isz-1;i++)
	{
		for(j=i+1;j<isz;j++)
		{
			if(v[i]<v[j])
				return false;
		}
	}
	return true;
}

template< class T > 
T Sum(std::vector<T>& v)
{
	int i;
	T val = T(0);
	for(i=0;i<v.size();i++) val += v[i];
	return val;
}

template< class T > 
T Sum(std::vector<T>& v,int iStart,int iEnd)
{
	int i;
	T val = T(0);
	for(i=iStart;i<v.size() && i<iEnd;i++) val += v[i];
	return val;
}

template< class T > 
T Avg(std::vector<T>& v)
{
	if(!v.size()) return T(0);
	int i;
	T val = T(0);
	for(i=0;i<v.size();i++) val += v[i];
	return val / (T) v.size();
}

template< class T > 
T Avg(std::vector<T>& v,int iStart,int iEnd)
{
	if(!v.size()) return T(0);
	int i;
	T val = T(0);
	for(i=iStart;i<iEnd;i++) val += v[i];
	return val / (T) (iEnd-iStart);
}

template< class T >
T Variance(std::vector<T>& v,int iStart,int iEnd)
{	T avg = Avg(v,iStart,iEnd);
	int i;
	T var = T(0) , val = T(0);
	for(i=iStart;i<iEnd;i++)
	{	val = v[i]-avg;
		var += val*val;
	}
	var /= (T)(iEnd-iStart);
	return var;
}

//returns sum(1..N)
inline int IntegerSum(int N)
{
	return (N*N+N)/2;
}

using namespace std;

//#define SAMPLE_MEAN

// Variance-covariance matrix: creation
template<class REAL>
void CovarMat(const vector< vector<REAL> >& vDataIn,int iRows,int iCols,vector< vector<REAL> >& vCovarMat, vector<REAL>& vMean,bool bCorrelMat = false)
{	// Create iCols * iCols covariance matrix from given iRow * iCol data matrix. 
	int i, j, x2, x, y, j1 , j2;
	// Allocate storage for mean vector 
	vMean = vector<REAL>(iCols,0.0);

	vector< vector<REAL> > vData(vDataIn); // copy it!

#ifdef SAMPLE_MEAN
	REAL N = iRows - 1;
#else
	REAL N = iRows;
#endif

	// Determine mean of row vectors of input data matrix 
	for(x=0;x<iCols;x++)
	{	for(y=0;y<iRows;y++)
			vMean[x] += vData[y][x];
		vMean[x] /= N;
	}

	vector<REAL> vStdev(iCols);
	
	// Center the row vectors.
	for(y=0;y<iRows;y++)
    {	for(x=0;x<iCols;x++)
		{	vData[y][x] -= vMean[x];
		}
    }

	// compute std-dev
	//const bool bCorrelMat = false;
	if(bCorrelMat) for(x=0;x<iCols;x++)
	{	for(y=0;y<iRows;y++)
			vStdev[x]+=vData[y][x]*vData[y][x];
		vStdev[x] /= N;
		vStdev[x] = sqrt(vStdev[x]);
	}

	vCovarMat = vector<vector<REAL> >(iCols,vector<REAL>(iCols,0.0));

	// Calculate the iCols * iCols covariance matrix.
	for(j1=0;j1<iCols;j1++)
	{	for(j2=0;j2<=j1;j2++)
		{	vCovarMat[j1][j2]=0.0;
			for(i=0;i<iRows;i++)
				vCovarMat[j1][j2] += ( vData[i][j1] * vData[i][j2] );
			if(bCorrelMat)	//correlation matrix
				vCovarMat[j1][j2] /= (N * vStdev[j1] * vStdev[j2] );
			else	//covariance matrix
				vCovarMat[j1][j2] /= N;
			vCovarMat[j2][j1]=vCovarMat[j1][j2];
		}
	}
}

// Variance-covariance matrix: creation
template<class REALD , class REALC, class REALM>
void CovarMat(  A2D<REALD>& vData , int iRows, int iCols, A2D<REALC>& vCovarMat, vector<REALM>& vMean,bool bCorrelMat = false)
{	// Create iCols * iCols covariance matrix from given iRow * iCol data matrix. 
	int i, x, y, j1 , j2;
	// Allocate storage for mean vector 
	vMean = vector<REALM>(iCols);

#ifdef SAMPLE_MEAN
	REALM N = iRows - 1;
#else
	REALM N = iRows;
#endif

	// Determine mean of row vectors of input data matrix 
	for(x=0;x<iCols;x++)
	{	vMean[x] = 0.0;
		for(y=0;y<iRows;y++)
		{	vMean[x] += vData[y][x];
		}
		vMean[x] /= N;
	}

	vector<REALM> vStdev(iCols);
	
	// Center the row vectors.
	for(y=0;y<iRows;y++)
    {	for(x=0;x<iCols;x++)
		{	vData[y][x] -= vMean[x];
		}
    }

	// compute std-dev
	//const bool bCorrelMat = false;
	if(bCorrelMat) for(x=0;x<iCols;x++)
	{	for(y=0;y<iRows;y++)
			vStdev[x]+=vData[y][x]*vData[y][x];
		vStdev[x] /= N;
		vStdev[x] = sqrt(vStdev[x]);
	}

	vCovarMat.Init(iCols,iCols);

	if(bCorrelMat)	// Calculate the iCols * iCols correlation matrix.
	{	for(j1=0;j1<iCols;j1++)
		{	for(j2=0;j2<=j1;j2++)
			{	vCovarMat[j1][j2]=0.0;
				for(i=0;i<iRows;i++)
				{	vCovarMat[j1][j2] += (vData[i][j1] * vData[i][j2]);
				}
				vCovarMat[j1][j2] /= ((REALC) N * vStdev[j1] * vStdev[j2] );
				vCovarMat[j2][j1]=vCovarMat[j1][j2];
			}
		}
	}
	else	// Calculate the iCols * iCols covariance matrix.
	{	for(j1=0;j1<iCols;j1++)
		{	for(j2=0;j2<=j1;j2++)
			{	vCovarMat[j1][j2]=0.0;
				for(i=0;i<iRows;i++)
				{	vCovarMat[j1][j2] += (vData[i][j1] * vData[i][j2]);
				}
				vCovarMat[j1][j2] /= N;
				vCovarMat[j2][j1]=vCovarMat[j1][j2];
			}
		}
	}

	// Un-Center the row vectors.
	for(y=0;y<iRows;y++)
    {	for(x=0;x<iCols;x++)
		{	vData[y][x] += vMean[x];
		}
    }
}

// Variance-covariance matrix: creation
template<class REAL>
void CovarMat(const vector< REAL >& vDataIn,int iRows,int iCols,vector< vector<REAL> >& vCovarMat, vector<REAL>& vMean,bool bCorrelMat = false)
{	// Create iCols * iCols covariance matrix from given iRow * iCol data matrix. 
	int i, x, y, j1 , j2;
	// Allocate storage for mean vector 
	vMean = vector<REAL>(iCols);

	vector< REAL > vData(vDataIn); // copy it!

#ifdef SAMPLE_MEAN
	REAL N = iRows - 1;
#else
	REAL N = iRows;
#endif

	// Determine mean of row vectors of input data matrix 
	for(x=0;x<iCols;x++)
	{	vMean[x] = 0.0;
		for(y=0;y<iRows;y++)
		{	vMean[x] += vData[y*iCols+x];
		}
		vMean[x] /= N;
	}

	vector<REAL> vStdev(iCols);
	
	// Center the row vectors.
	for(y=0;y<iRows;y++)
    {	for(x=0;x<iCols;x++)
		{	vData[y*iCols+x] -= vMean[x];
		}
    }

	// compute std-dev
	//const bool bCorrelMat = false;
	if(bCorrelMat) for(x=0;x<iCols;x++)
	{	for(y=0;y<iRows;y++)
			vStdev[x]+=vData[y*iCols+x]*vData[y*iCols+x];
		vStdev[x] /= N;
		vStdev[x] = sqrt(vStdev[x]);
	}

	vCovarMat = vector<vector<REAL> >(iCols,vector<REAL>(iCols,0.0));

	// Calculate the iCols * iCols covariance matrix.
	for(j1=0;j1<iCols;j1++)
	{	for(j2=0;j2<=j1;j2++)
		{	vCovarMat[j1][j2]=0.0;
			for(i=0;i<iRows;i++)
			{	vCovarMat[j1][j2] += (vData[i*iCols+j1] * vData[i*iCols+j2]);
			}
			if(bCorrelMat)	//correlation matrix
				vCovarMat[j1][j2] /= ( N * vStdev[j1] * vStdev[j2] );
			else	//covariance matrix
				vCovarMat[j1][j2] /= N;
			vCovarMat[j2][j1]=vCovarMat[j1][j2];
		}
	}
}

template< class T >
int GetMaxIndex(vector<T>& v)
{
	int iSz = v.size(), i = 0, iMaxIDX = 0;
	if(iSz<1) return -1;
	T maxVal = v[0];
	for(i=1;i<iSz;i++)
	{
		if(v[i]>maxVal)
		{
			maxVal=v[i];
			iMaxIDX=i;
		}
	}
	return iMaxIDX;
}

template< class T >
int GetMinIndex(vector<T>& v)
{
	int iSz = v.size(), i = 0, iMinIDX = 0;
	if(iSz<1) return -1;
	T minVal = v[0];
	for(i=1;i<iSz;i++)
	{
		if(v[i]<minVal)
		{
			minVal=v[i];
			iMinIDX=i;
		}
	}
	return iMinIDX;
}

template< class T >
T Energy(vector<T>& v)
{
	int iSz = v.size() , i = 0;
	T sum(0);
	if(iSz<1) return sum;
	for(;i<iSz;i++) sum += v[i]*v[i];
	sum = sqrt(sum);
	return sum / (double) iSz;
}

//chi-squared cumulative distribution function
//returns probability of chi-squared with df degrees of freedom has value <= x
double chisqCDF(double x, double df);

inline bool PowOf2(int N)
{
	if(N < 0) return false;
	while(N>1)
	{
		if(N%2==1) return false;
		N /= 2;
	}
	return true;
}

inline int CeilPowOf2(int N)
{
	int T = 1;
	while(T < N) T *= 2;
	return T;
}

