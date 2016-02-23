// $Id: Vertex.h,v 1.3 2011/01/08 01:14:24 samn Exp $ 
// Vertex.h: interface for the CVertex class.
//
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
//////////////////////////////////////////////////////////////////////

#ifndef VERTEX_H
#define VERTEX_H

#include "MyObj.h"
#include "A2D.h"
#include "Cluster.h"
#include "Log.h"
#include <deque>
#include <map>

using namespace std;

////////////////////////////////////////////////////////////////////////
// CVertex
class CVertex : public CMyObject  
{
	friend class CVerxStack;
protected:
	VERTEX	m_Vertex;		// Main values of vector
	char	m_Clust;//original cluster
public:
	CVertex(){m_Clust=0;};
	virtual ~CVertex(){};
	int		GetClust(){ return m_Clust; }
	void    SetClust(char c) { m_Clust = c; }
	//gets number of dims in this vector, remember that dim 0 is the # of clusters this vector belongs to
	int		GetNumDims() { return m_Vertex.size(); };
	float	GetValue(int mIndex){ return *(m_Vertex.begin()+mIndex);};
	void	SetValue(int mIndex,float Value){ m_Vertex[mIndex]=Value; };
	void AddPnt(float);
	friend class CCluster;
};

double* getrank (int n, double data[]);

////////////////////////////////////////////////////////////////////////
// CVerxStack
class CVerxStack : public CMyObject
{
public:
	// stacks
	MY_STACK		m_VerxStack;	// main vectors of spikes	
	int				m_iNumDim;
	int				m_NumVerx;		// in memory
	VERTEX			m_MainMin,m_MainMax; // without noise
	VERTEX			m_MainNoisedMin,m_MainNoisedMax;
	VERTEX			m_MainRange,m_MainStdev,m_MainMean,m_MainEntropy;
	bool			m_bNormFloatV;//whether to normalize dimensions to be between 0-1 when computing KLD
	CCluster		m_oClusters;

	inline char GetVClust(CVertex* verx)
	{
    	  return verx->GetClust();
	}

	//get number of non-noise vertices in each cluster
	//uses whichDraw to determine which mode to get counts for
	void GetCounts(vector<int>& vCounts,int iClusts)
	{
		vCounts = vector<int>(iClusts+1);
		MY_STACK::iterator IT = m_VerxStack.begin();
		for(;IT!=m_VerxStack.end();IT++)
		{
			CVertex* verx = (CVertex*)*IT;
			vCounts[GetVClust(verx)]++;
		}
	}

	//get clust IDs
	void GetClustIDs(vector<int>& vIDs)
	{
		vIDs = vector<int>(NumNonNoiseVertices());
		int iV = 0;
		MY_STACK::iterator IT = m_VerxStack.begin();
		for(;IT!=m_VerxStack.end();IT++)
		{
			CVertex* verx = (CVertex*)*IT;
			vIDs[iV++] = verx->GetClust();
		}
	}

	inline int NumNonNoiseVertices()
	{
		return m_VerxStack.size();
	}

	//get 2D vector of vertex values NOT normalized between 0 - 1
	template < class T > 
	T** GetPV(int& iRows,int& iCols)
	{
		MY_STACK::iterator Index;
		int iSz = m_VerxStack.size();
		int iDims = m_iNumDim;
		T** p = Allocate2DArray<T>(iSz,iDims);
		iRows = iSz; iCols = iDims;
		int iV=0;
		for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++)
		  {
		    CVertex* verx = (CVertex*) *Index;
		    int iD = 0;
		    for(iD=0;iD<iDims;iD++)
		      {
			p[iV][iD] = verx->GetValue(iD);
		      }
		    iV++;
		  }
		return p;
	}

	//get 2D vector of vertex values as indices to bins in a distrib
	bool GetVertexFloatps(int& iRows,int& iCols,int iBins, vector< vector<float>* >& vFloatps)
	{
		MY_STACK::iterator Index;
		int iSz = m_VerxStack.size();		
		int iDims = m_iNumDim;
		vFloatps = vector< vector<float>* >(iSz);
		int iV = 0;
		for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++)
		  {
		    CVertex* verx = (CVertex*) *Index;
		    vFloatps[iV] = &verx->m_Vertex;
		  }		
		iRows = iSz; iCols = iDims;
		return true;
	}

	bool GetFloatV(int& iRows,int& iCols,vector<float>& vFloat, vector<float>& vRange)
	{	
		MY_STACK::iterator Index;
		int iSz = NumNonNoiseVertices();
		
		//just gets # of dimensions without x,y,time locations
		int iDims = m_iNumDim;

		CalcDimStats();
		iRows = iSz;
		iCols = iDims;

		vFloat = vector<float>(iSz*iCols);

		Write2Log("m_bNormFloatV = %s",m_bNormFloatV?"true":"false");

		int i = 0;
		vRange = vector<float>(iDims);
		for(i=0;i<iDims;i++) vRange[i]=GetMax(i)-GetMin(i);

		if(m_bNormFloatV)	//do normalization of data
		{	int iV=0, j = 0, iD = 0;
			for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++)
			{	CVertex* verx = (CVertex*) *Index;
				for(iD=0;iD<iDims;iD++)					
					vFloat[j++]=(verx->GetValue(iD) - GetMin(iD)) / vRange[iD];
			}
		}
		else	//don't do normalization
		{	int iV=0, j = 0, iD = 0;
			for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++)
			{	CVertex* verx = (CVertex*) *Index;
				for(iD=0;iD<iDims;iD++)				
					vFloat[j++]=verx->GetValue(iD);
			}
		}
		return true;
	}

	bool GetFloat2D(int& iRows,int& iCols,A2D<float>& vFloat, vector<float>& vRange)
	{	
		MY_STACK::iterator Index;
		int iSz = NumNonNoiseVertices();
		
		//just gets # of dimensions without x,y,time locations
		int iDims = m_iNumDim;

		CalcDimStats();
		iRows = iSz;
		iCols = iDims;

		vFloat.Init(iSz,iCols);

		Write2Log("m_bNormFloatV = %s",m_bNormFloatV?"true":"false");

		int i = 0, Y = 0;
		vRange = vector<float>(iDims+1);
		for(i=0;i<iDims;i++) vRange[i]=GetMax(i)-GetMin(i);

		if(m_bNormFloatV)	//do normalization of data
		{	int iV=0, j = 0, iD = 0;
			for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++)
			{	CVertex* verx = (CVertex*) *Index;
				j=0;
				for(iD=0;iD<iDims;iD++)					
					vFloat[iV][j++]=(verx->GetValue(iD) - GetMin(iD)) / vRange[iD];
			}
		}
		else	//don't do normalization
		{	int iV=0, j = 0, iD = 0;
			for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++)
			{	CVertex* verx = (CVertex*) *Index;
				j=0;
				for(iD=0;iD<iDims;iD++)				
					vFloat[Y][j++]=verx->GetValue(iD);
			}
		}
		return true;
	}


	//get 2D vector of vertex values as indices to bins in a distrib
	bool GetVBinIDs(int& iRows,int& iCols,int iBins, vector< vector<int> >& vBinIDs)
	{
		MY_STACK::iterator Index;
		int iSz = m_VerxStack.size();		
		int iDims = m_iNumDim;
		vBinIDs = vector< vector<int> >(iSz, vector<int>(iDims));
		iRows = iSz; iCols = iDims;
		int iV=0, iB=0;
		for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++)
		{
			CVertex* verx = (CVertex*) *Index;
			int iD = 0;
			for(iD=0;iD<iDims;iD++)
			{	//+1 since index 0 is # of clusters vertex is in
				iB = (iBins)*(verx->GetValue(iD) - GetMin(iD)) / (GetMax(iD) - GetMin(iD));
				if(iB >= iBins) iB = iBins - 1;
				vBinIDs[iV][iD] = iB;
			}
		}
		return true;
	}

	//get 2D vector of vertex values as indices to bins in a distrib
	//if bExcludeNoise, no noise vertices will be stored and user
	//should know the indices will be offset accordingly
	template < class T > 
	T** GetVBinIDs(int& iRows,int& iCols,int iBins,bool bRank=false)
	{
		MY_STACK::iterator Index;
		int iSz = m_VerxStack.size();
		int iDims = m_iNumDim;
		T** p = Allocate2DArray<T>(iSz,iDims);
		iRows = iSz; iCols = iDims;
		int iV=0,iB=0,iNNSz=0;
		if(bRank)
		{
			iNNSz = NumNonNoiseVertices();
			double* ptmp = (double*) malloc(sizeof(double)*iNNSz);
			int iD;
			for(iD=0;iD<iDims;iD++)
			{	iV=0;
				for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++)
				{
					CVertex* verx = (CVertex*) *Index;
					ptmp[iV] = verx->GetValue(iD);
					iV++;
				}
				double* prank = getrank(iNNSz,ptmp);
				for(iV=0;iV<iNNSz;iV++) 
				{	iB = iBins * prank[iV] / (iNNSz-1);
					if(iB >= iBins) iB = iBins - 1;
					p[iV][iD] = iB;
				}
				free(prank);
			}
			free(ptmp);
		}
		else
		{
			for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++)
			{
				CVertex* verx = (CVertex*) *Index;
				int iD = 0;
				for(iD=0;iD<iDims;iD++)
				{
					iB = (iBins)*(verx->GetValue(iD) - GetMin(iD)) / (GetMax(iD) - GetMin(iD));
					if(iB >= iBins) iB = iBins - 1;
					p[iV][iD] = iB;
				}
				iV++;
			}
		}
		return p;
	}

	//get 2D vector of vertex values NOT normalized between 0 - 1
	template < class T > 
	void GetV( std::vector<std::vector<T> >& V)
	{
		MY_STACK::iterator Index;
		int iSz = m_VerxStack.size();
		int iDims = m_iNumDim;
		V = vector< vector<T> >(iSz, vector<T>(iDims));
		int iV=0;
		for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++)
		{
			CVertex* verx = (CVertex*) *Index;
			int iD = 0;
			for(iD=0;iD<iDims;iD++)
				V[iV][iD] = verx->GetValue(iD);
			iV++;
		}
	}

	//get 2D vector of vertex values normalized between 0 - 1
	template < class T >
	T** NormalizedPV(int& iRows,int& iCols)
	{
		MY_STACK::iterator Index;
		int iSz = m_VerxStack.size();
		int iDims = m_iNumDim;
		T** p = Allocate2DArray<T>(iSz,iDims);
		iRows = iSz; iCols = iDims;
		vector<T> vRange(iDims);
		int i;
		for(i=0;i<iDims;i++) vRange[i] = GetMax(i) - GetMin(i);
		int iV=0;
		for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++)
		{
			CVertex* verx = (CVertex*) *Index;
			int iD = 0;
			for(iD=0;iD<iDims;iD++)
				p[iV][iD] = (verx->GetValue(iD) - GetMin(iD)) /  vRange[iD];;
			iV++;
		}
		return p;
	}

	//get 2D vector of vertex values normalized between 0 - 1
	template < class T >
	void NormalizedV(std::vector< std::vector<T> >& vN)
	{
		MY_STACK::iterator Index;
		int iSz = m_VerxStack.size();
		int iDims = m_iNumDim;
		vN = vector< vector<T> >(iSz, vector<T>(iDims));
		int i;
		vector<T> vRange(iDims);
		for(i=0;i<iDims;i++) vRange[i] = GetMax(i) - GetMin(i);
		int iV=0;
		for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++)
		{
			CVertex* verx = (CVertex*) *Index;
			int iD = 0;
			for(iD=0;iD<iDims;iD++)
				vN[iV][iD] = (verx->GetValue(iD) - GetMin(iD)) /  vRange[iD];
			iV++;
		}
	}
// Methods
protected:
	void	CalcAfterLoad();
public:		
	CVerxStack();
	~CVerxStack(){ SetEmpty(); };
	void	AddVrx(CMyObject *toStore);
	void	CalcDimStats();	//mean, stdev, min, max, range, entropy
	void	CalcMinMax();	
	int		GetClust(int NoOfVerx);
	int		GetDimension(){ return m_iNumDim; };
	int		GetNumVerx(){ return m_NumVerx; };
	float	GetMin(int Index);
	float	GetMax(int Index);
	void	SetEmpty();	
	bool ReadAsciiData(char* path);
	friend class CCluster;
private:
	void InitClust(int iMaxClust);
};

#endif
