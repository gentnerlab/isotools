// $Id: IsolationDist.cpp,v 1.5 2008/08/11 20:55:44 samn Exp $ 

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


#include "IsolationDist.h"
#include "Cluster.h"
#include "nr.h"
#include "Log.h"
#include <algorithm>
#include "WCMath.h"

using namespace std;

typedef double ID_T; //need double for extra precision in matrix operations

void MatMult(A2D<ID_T>& a, A2D<ID_T>& b,A2D<ID_T>& out)
{	int iRowsA = 0, iColsA = 0, iRowsB = 0 , iColsB = 0;
	try
	{
		iRowsA = a.Rows(); if(!iRowsA) return; iColsA = a.Cols();
		iRowsB = b.Rows(); if(!iRowsB) return; iColsB = b.Cols();
		if(iColsA!=iRowsB)
			return;
		out.Init(iRowsA,iColsB); out.Fill(0.0);
		int i,j,k;
		for(i=0;i<iRowsA;i++)
			for(j=0;j<iColsB;j++)
				for(k=0;k<iColsA;k++)
					out[i][j] += a[i][k] * b[k][j];
	}
	catch(...)
	{
		Write2Log("Exception in MatMult iRowsA=%d iColsA=%d iRowsB=%d iColsB=%d",iRowsA,iColsA,iRowsB,iColsB);
	}
}

void MatMult(vector<vector<ID_T> >& a, vector<vector<ID_T> >& b,vector<vector<ID_T> >& out)
{	int iRowsA = 0, iColsA = 0, iRowsB = 0 , iColsB = 0;
	try
	{
		iRowsA = a.size(); if(!iRowsA) return; iColsA = a[0].size();
		iRowsB = b.size(); if(!iRowsB) return; iColsB = b[0].size();
		if(iColsA!=iRowsB)
			return;
		out = vector<vector<ID_T> >(iRowsA,vector<ID_T>(iColsB));
		int i,j,k;
		for(i=0;i<iRowsA;i++)
			for(j=0;j<iColsB;j++)
				for(k=0;k<iColsA;k++)
					out[i][j] += a[i][k] * b[k][j];
	}
	catch(...)
	{
		Write2Log("Exception in MatMult iRowsA=%d iColsA=%d iRowsB=%d iColsB=%d",iRowsA,iColsA,iRowsB,iColsB);
	}
}

bool Row2Col(vector<vector<float> >& vDataIn,int iRowID,vector<vector<float> >& vDataOut,int iColID)
{
	if(!vDataIn.size())
		return false;
	if(!vDataOut.size())
		return false;
	int iCols = vDataIn[0].size();
	int x;
	for(x=0;x<iCols;x++)
	{
		vDataOut[x][iColID]=vDataIn[iRowID][x];
	}
	return true;
}

bool Transpose(vector<vector<float> >& vDataIn,vector<vector<float> >& vDataOut)
{
	int iRows = vDataIn.size();
	int iCols = vDataIn[0].size();
	if(!iRows || !iCols) return false;
	vDataOut = vector<vector<float> >(iCols,vector<float>(iRows));
	int x,y;
	for(y=0;y<iRows;y++)
	{	for(x=0;x<iCols;x++)
		{
			vDataOut[x][y]=vDataIn[y][x];
		}
	}
	return true;
}

bool InvertMatrix(vector< vector<ID_T> >& vData)
{
	ID_T d;
	vector<vector<ID_T> > vTmp(vData);
	int i,j;
	int N = vData.size();
	vector<int> indx(N);
	vector<ID_T> col(vData[0].size());
	if(!NR::ludcmp(vTmp,indx,d))//Decompose the matrix just once.
	{	Write2Log("Couldn't invert matrix!");
		return false;		
	}
	for(j=0;j<N;j++) 
	{ //Find inverse by columns.
		for(i=0;i<N;i++)
			col[i]=0.0;
		col[j]=1.0;
		NR::lubksb(vTmp,indx,col);
		for(i=0;i<N;i++) 
			vData[i][j]=col[i];
	}

	return true;
}

bool InvertMatrix(A2D<ID_T>& vData)
{
	ID_T d;
	A2D<ID_T> vTmp(vData.Rows(),vData.Cols());
	int i,j;
	for(i=0;i<vData.Rows();i++)
		for(j=0;j<vData.Cols();j++)
			vTmp[i][j]=vData[i][j];
	
	int N = vData.Rows();
	vector<int> indx(N);
	vector<ID_T> col(vData.Cols());
	if(!NR::ludcmp(vTmp,indx,d))//Decompose the matrix just once.
	{	Write2Log("Couldn't invert matrix!");
		return false;		
	}
	for(j=0;j<N;j++) 
	{ //Find inverse by columns.
		for(i=0;i<N;i++)
			col[i]=0.0;
		col[j]=1.0;
		NR::lubksb(vTmp,indx,&col[0]);
		for(i=0;i<N;i++) 
			vData[i][j]=col[i];
	}

	return true;
}

// Variance-covariance matrix: creation
/*
void CovarMat(const vector< vector<ID_T> >& vDataIn,int iRows,int iCols,vector< vector<ID_T> >& vCovarMat, vector<ID_T>& vMean)
{	// Create iCols * iCols covariance matrix from given iRow * iCol data matrix. 
	int i, j, x2, x, y, j1 , j2;
	// Allocate storage for mean vector 
	vMean = vector<ID_T>(iCols);

	vector< vector<ID_T> > vData = vDataIn; // copy it!

	// Determine mean of row vectors of input data matrix 
	for(x=0;x<iCols;x++)
	{	vMean[x] = 0.0;
		for(y=0;y<iRows;y++)
		{	vMean[x] += vData[y][x];
		}
		vMean[x] /= (ID_T)iRows;
	}

	vector<ID_T> vStdev(iCols);
	
	// Center the row vectors.
	for(y=0;y<iRows;y++)
    {	for(x=0;x<iCols;x++)
		{	vData[y][x] -= vMean[x];
		}
    }

	// compute std-dev
	const bool bCorrelMat = false;
	if(bCorrelMat) for(x=0;x<iCols;x++)
	{	for(y=0;y<iRows;y++)
			vStdev[x]+=vData[y][x]*vData[y][x];
		vStdev[x] /= (ID_T) iRows;
		vStdev[x] = sqrt(vStdev[x]);
	}

	vCovarMat = vector<vector<ID_T> >(iCols,vector<ID_T>(iCols,0.0));

	// Calculate the iCols * iCols covariance matrix.
	for(j1=0;j1<iCols;j1++)
	{	for(j2=0;j2<=j1;j2++)
		{	vCovarMat[j1][j2]=0.0;
			for(i=0;i<iRows;i++)
			{	vCovarMat[j1][j2] += (vData[i][j1] * vData[i][j2]);
			}
			if(bCorrelMat)	//correlation matrix
				vCovarMat[j1][j2] /= ((ID_T) iRows * vStdev[j1] * vStdev[j2] );
			else	//covariance matrix
				vCovarMat[j1][j2] /= (ID_T) (iRows); // -1);
			vCovarMat[j2][j1]=vCovarMat[j1][j2];
		}
	}
}
*/
//gets the inverse of the covariance matrix for a cluster
//input vFloat has each row as a data vector, each column as a dimension, to get data vector i do vFloat[i*iCols]
//bool Clust2CovarMatInv(vector<vector< ID_T > >& vCovarMat,vector<ID_T>& vMean, vector<int>& vClustIDs,int iClustID,vector<double>& vFloat,int iRows,int iCols,int& iClustSz)
bool Clust2CovarMatInv(A2D< ID_T >& vCovarMat,vector<ID_T>& vMean, vector<int>& vClustIDs,int iClustID,vector<double>& vFloat,int iRows,int iCols,int& iClustSz)
{
	int i = 0, j = 0;
		iClustSz = count(vClustIDs.begin(),vClustIDs.end(),iClustID);
	if(!iClustSz)
		return false;
	A2D<double> vClustData(iClustSz,iCols);
	vClustData.Fill(0.0);
	for(i=0;i<iRows;i++)
		if(vClustIDs[i]==iClustID)
			copy(&vFloat[i*iCols],&vFloat[i*iCols+iCols],vClustData[j++]);

	Write2LogPlain("clust %d sz=%d\n",iClustID,iClustSz);
	
	CovarMat(vClustData,iClustSz,iCols,vCovarMat,vMean);
		
	if(false){
	Write2LogPlain("clust %d covar mat %dX%d\n",iClustID,vCovarMat.Rows(),vCovarMat.Cols());
	for(i=0;i<vCovarMat.Rows();i++){
		for(j=0;j<vCovarMat.Cols();j++){
			Write2LogPlain("%g ",vCovarMat[i][j]);
		}
		Write2LogPlain("\n");
	}
	}

#ifdef _DEBUG
	A2D<ID_T> mtmp(vCovarMat);
	
	bool b = InvertMatrix(vCovarMat);

	Write2LogPlain("inv covar mat");
	for(i=0;i<vCovarMat.Rows();i++){
		for(j=0;j<vCovarMat.Cols();j++){
			Write2LogPlain("%g ",vCovarMat[i][j]);
		}
		Write2LogPlain("\n");
	}
	
	A2D<ID_T> vident;
	MatMult(mtmp,vCovarMat,vident);
	
	Write2LogPlain("is this identity matrix?");
	for(i=0;i<vident.Rows();i++){
		for(j=0;j<vident.Cols();j++){
			Write2LogPlain("%g ",vident[i][j]);
		}
		Write2LogPlain("\n");
	}
	
	return b;
#else
	return InvertMatrix(vCovarMat);
#endif
}

ID_T DotProd(vector<ID_T>& v1,vector<ID_T>& v2)
{
	int iSz = v1.size(), i = 0;
	ID_T val = 0.0;
	for(i=0;i<iSz;i++) val += v1[i]*v2[i];
	return val;
}

//squared mahalanobis distance of pData to vMean using inverse of covariance matrix, vCovarMatInv
ID_T MahalDistSQ(double* pData,int iRowSz,vector< vector<ID_T> >& vCovarMatInv, vector<ID_T>& vMean)
{	// (x-mean)'*covar^-1*(x-mean)
	vector<ID_T> vRow(iRowSz,0.0), vTmp(iRowSz,0.0);
	int i, iCols = iRowSz , j = 0;
	for(i=0;i<iCols;i++)	//center row data vector
		vRow[i]=pData[i] - vMean[i];
	
	//multiply vector by matrix
	for(i=0;i<iCols;i++)
		for(j=0;j<iCols;j++)
			vTmp[i] += vRow[j] * vCovarMatInv[j][i];

	//dot product
	return DotProd(vTmp,vRow);
}

ID_T MahalDistSQ(double* pData,int iRowSz,A2D<ID_T>& vCovarMatInv, vector<ID_T>& vMean)
{	// (x-mean)'*covar^-1*(x-mean)
	vector<ID_T> vRow(iRowSz,0.0), vTmp(iRowSz,0.0);
	int i, iCols = iRowSz , j = 0;
	for(i=0;i<iCols;i++)	//center row data vector
		vRow[i]=pData[i] - vMean[i];
	
	//multiply vector by matrix
	for(i=0;i<iCols;i++)
		for(j=0;j<iCols;j++)
			vTmp[i] += vRow[j] * vCovarMatInv[j][i];

	//dot product
	return DotProd(vTmp,vRow);
}

float IsolationDist(vector<int>& vClustIDs,int iClustID,vector<double>& vFloat,int iRows,int iCols)
{
	try
	{
		A2D<double> vCovarMatInv;
		vector<ID_T> vMean;
		int iClustSz = 0;
		if(!Clust2CovarMatInv(vCovarMatInv,vMean,vClustIDs,iClustID,vFloat,iRows,iCols,iClustSz))
			return -1.0;
		if(iClustSz*2>=iRows || iClustSz<2) //need at least as many points outside of cluster as within it
			return -1.0;

		//Write2Log("vMean for cluster %d : ",iClustID); WriteVec2Log(vMean);

		vector<ID_T> vDists; 
		int i;
		for(i=0;i<iRows;i++)
		{	if(vClustIDs[i]==iClustID)
				continue;
			ID_T dDist = MahalDistSQ(&vFloat[i*iCols],iCols,vCovarMatInv,vMean);
			if(dDist>=0.0)
				vDists.push_back(dDist);
			else
				Write2Log("Negative MahalDistSq = %.4f",dDist);
		}

		sort(vDists.begin(),vDists.end());

		if(false){
		Write2LogPlain("clust %d mahal distances:\n",iClustID);
		for(i=0;i<vDists.size();i++){ Write2LogPlain("%g ",vDists[i]); if(i%20==0) Write2LogPlain("\n"); }
		Write2LogPlain("\n");}

	#ifdef _DEBUG
		Write2Log("meandist=%.4f mediandist=%.4f",Avg(vDists),vDists[vDists.size()/2]);
	#endif

		if(iClustSz-1>=0 && iClustSz-1<vDists.size())
			return vDists[iClustSz-1];
		return -1.0;
	}
	catch(...)
	{
		Write2Log("Exception in IsolationDist!");
		return -1.0;
	}
}

float LRatio(vector<int>& vClustIDs,int iClustID,vector<double>& vFloat,int iRows,int iCols)
{
	try
	{
		A2D<ID_T> vCovarMatInv;
		vector<ID_T> vMean;
		int iClustSz = 0;
		if(!Clust2CovarMatInv(vCovarMatInv,vMean,vClustIDs,iClustID,vFloat,iRows,iCols,iClustSz))
			return -1.0;
		if(iClustSz<1 || iRows-iClustSz<1)
			return -1.0;

		int i;
		ID_T L = 0.0;
		for(i=0;i<iRows;i++)
		{	if(vClustIDs[i]==iClustID)
				continue;
			ID_T dDist = MahalDistSQ(&vFloat[i*iCols],iCols,vCovarMatInv,vMean);
			if(dDist>=0.0)
				L += (1.0 - chisqCDF(dDist,iCols));
			else
				Write2Log("Negative MahalDistSq = %.4f",dDist);
		}
		Write2Log("Clust%d : L = %g , L-Ratio = %g", iClustID, L , L / iClustSz);
		return L / iClustSz;
	}
	catch(...)
	{
		Write2Log("Exception in LRatio!");
		return -1.0;
	}
}
