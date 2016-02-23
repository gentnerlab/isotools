// $Id: Hist.h,v 1.2 2011/01/08 01:13:48 samn Exp $ 
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

#ifndef HIST_H_
#define HIST_H_

#include <map>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "Cluster.h"
#include "WCMath.h"
#include "KDTree.h"
#include "Log.h"

//#define DO_TIMING

//binary logarithm
inline prob_t log2(prob_t d)
{
	static prob_t dl2 = log(2.0);
	return log(d) / dl2;
}

inline double log2(double d)
{
	static double dl2 = log(2.0);
	return log(d) / dl2;
}

prob_t Prob(int,int);

void InitProbs(int iMaxNumElems);

struct ProbInitFree
{
	ProbInitFree(int i);
	~ProbInitFree();
};

inline bool GZeroMinElem(float f1,float f2)
{
	return f1 > 0.0 && (f1 < f2 || f2 <= 0.0);
}

inline float GZeroMinElem(vector<float>& v)
{
	float fm = *std::max_element(v.begin(),v.end());
	int i=0,sz=v.size();
	for(;i<sz;i++)
		if(v[i]<fm && v[i]>0.0)
			fm=v[i];
	return fm;
}

inline float MinElem(vector<float>& v,bool bAllowZero)
{
	float fm = *std::max_element(v.begin(),v.end());
	int i=0,sz=v.size();
	for(;i<sz;i++)
		if(v[i]<fm && (bAllowZero || v[i]>0.0))
			fm=v[i];
	return fm;
}

inline int MinIdx(vector<float>& v,bool bAllowZero)
{
	float fm = *std::max_element(v.begin(),v.end());
	int i=0,sz=v.size(), minID = 0;
	bool bFound = false;
	for(;i<sz;i++)
	{
		if(v[i]<fm && (bAllowZero || v[i]>0.0))
		{
			bFound = true;
			fm=v[i];
			minID=i;
		}
	}
	if(!bFound) return -1;
	return minID;
}

struct Neighbor
{
	prob_t m_dist;
	int m_id;
	Neighbor(int id,prob_t dist)
		:m_id(id),m_dist(dist){}
	Neighbor()
		:m_id(0),m_dist(0){}
};

inline bool operator<(Neighbor& n1,Neighbor& n2)
{
	return n1.m_dist<n2.m_dist;
}

inline bool operator<(const Neighbor& n1,const Neighbor& n2)
{
	return n1.m_dist<n2.m_dist;
}

//kd tree from biopython used for continuous probability distribution estimates 
class KDTreeHist
{
	NSKDTree::KDTree* m_pTree;

	int m_iDims;
	int m_iNumElems;	
	prob_t m_dPiPow;
	prob_t m_dTop;
	prob_t m_dGamma;

	vector<float> m_vData;
	
	vector<float> m_vProbs;

public:
	
	KDTreeHist()
		:m_iDims(0),
		 m_iNumElems(0),
		 m_pTree(0)
	{
	}

	virtual ~KDTreeHist()
	{
		if(m_pTree) delete m_pTree;
	}

	prob_t Top(){ return m_dTop; }
	prob_t PiPow(){ return m_dPiPow; }
	prob_t Gam(){ return m_dGamma; }

	bool SetData(vector<float>& vFloat,vector<int>& vIDs,int iNumPoints,int iDims,int iCID,bool bNot=false)
	{
		if(!vFloat.size() || iNumPoints<1) return false;
		int iV = 0;
		m_vData = vector<float>(iNumPoints*iDims);
		int j = 0;
		if(bNot)
		{
			for(iV=0;iV<vIDs.size();iV++)
			{
				if(vIDs[iV] != iCID && vIDs[iV]<1000)
				{
					int iD = 0;
					for(iD=0;iD<iDims;iD++)
					{
						m_vData[j++]=vFloat[iV*iDims+iD];
					}
				}
			}
		}
		else
		{
			for(iV=0;iV<vIDs.size();iV++)
			{
				if(vIDs[iV] == iCID)
				{
					int iD = 0;
					for(iD=0;iD<iDims;iD++)
					{
						m_vData[j++]=vFloat[iV*iDims+iD];
					}
				}
			}
		}
		m_iDims = iDims;		
		if(m_pTree) delete m_pTree;		
		m_pTree = new NSKDTree::KDTree(m_iDims,8,false);
		m_iNumElems = iNumPoints;
		const prob_t PI=3.14159265358979323846;
		m_dPiPow = pow((prob_t)PI,(prob_t)m_iDims/two);
		m_dTop = (one/(m_iNumElems-one));
		m_dGamma = Gamma(m_iDims/two+one);		
		m_pTree->set_data(&m_vData[0],iNumPoints);
		return true;
	}

	bool SetData(int iDims,float* pData,int iNumPoints)
	{
		if(!pData || iNumPoints<1) return false;

		m_iDims = iDims;
		if(m_pTree) delete m_pTree;
		m_pTree = new NSKDTree::KDTree(m_iDims,8,false);
		m_iNumElems = iNumPoints;

		m_vData = vector<float>(iNumPoints*m_iDims);
		memcpy(&m_vData[0],pData,iNumPoints*m_iDims*sizeof(float));

		const prob_t PI=3.14159265358979323846;

		m_dPiPow = pow((prob_t)PI,(prob_t)m_iDims/two);
		m_dTop = (one/(m_iNumElems-one));
		m_dGamma = Gamma(m_iDims/two+one);
		
		m_pTree->set_data(&m_vData[0],iNumPoints);
		
		return true;
	}

	int NumElems()
	{
		return m_iNumElems;
	}

	int NumDims()
	{
		return m_iDims;
	}

	float* operator[](int i)
	{
		if(i >= m_iNumElems) return 0;		
		return &m_vData[i*m_iDims];
	}

	float GetKNN(float* p,vector<Neighbor>& vnn,int iNNToFind,double fRadStart=8.0,double fRadFctr=2.0)
	{
		if(m_iNumElems == 1) return 0.0;
		
		float fRad = fRadStart;

		while(true)
		{
			m_pTree->search_center_radius_sq(p,fRad,iNNToFind);
			int iCount = m_pTree->get_count();
			if(iCount>iNNToFind)// || (iCount&&iNNToFind==1))
			{
				vector<float> vRadii(iCount);
				m_pTree->copy_radii_sq(&vRadii[0]);

				vector<long> vID(iCount);
				m_pTree->copy_indices(&vID[0]);

				vnn=vector<Neighbor>(iCount);

#ifdef _DEBUG
				Write2Log("found %d neighbors want %d",iCount,iNNToFind);
#endif

				int i = 0, j = 0;
				for(i=0;i<iCount;i++)
					if(vRadii[i]>0.0)
						vnn[j++] = Neighbor(vID[i],vRadii[i]);
				
				vnn.resize(j);
				std::sort(vnn.begin(),vnn.end());

				if(vnn.size()>iNNToFind)
				{
					vector<Neighbor> vNNTmp(iNNToFind);
					std::copy(vnn.begin(),vnn.begin()+iNNToFind,vNNTmp.begin());
					vnn=vNNTmp;
				}

#ifdef _DEBUG
				prob_t ttt = vnn[0].m_dist;
#endif
				
				return fRad;
			}
			//increase search radius
			fRad *= fRadFctr;
		}		
	}

	float GetKNN(float* p,Neighbor* vnn,int iNNToFind,double fRadStart,double fRadFctr)
	{
		if(m_iNumElems == 1) return 0.0;
		
		float fRad = fRadStart;

		while(true)
		{
			m_pTree->search_center_radius_sq(p,fRad,iNNToFind);
			int iCount = m_pTree->get_count();
			if(iCount>iNNToFind)// || (iCount&&iNNToFind==1))
			{
				vector<float> vRadii(iCount);
				m_pTree->copy_radii_sq(&vRadii[0]);

				vector<long> vID(iCount);
				m_pTree->copy_indices(&vID[0]);

				vector<Neighbor> vnnTMP(iCount);

				int i = 0, j = 0;
				for(i=0;i<iCount;i++)
					if(vRadii[i]>0.0)
						vnnTMP[j++] = Neighbor(vID[i],vRadii[i]);
				
				vnnTMP.resize(j);
				std::sort(vnnTMP.begin(),vnnTMP.end());

				for(i=0;i<iNNToFind;i++) vnn[i]=vnnTMP[i];

				return fRad;
			}
			//increase search radius
			fRad *= fRadFctr;
		}		
	}

	float GetAllKNN(A2D<Neighbor>& vnn,int iNNToFind,double fRadStart,double fRadFctr,vector<int>& vFound)
	{
		if(m_iNumElems == 1) return 0.0;
		
		float fRad = fRadStart;

		vnn.Fill(Neighbor(0,INF));
		vFound.resize(m_iNumElems);

		while(true)
		{
			m_pTree->neighbor_search_sq(fRadStart);
			int iCount = m_pTree->neighbor_get_count();
			if(iCount/2>m_iNumElems*iNNToFind)
			{
				Neighbor** pnn = vnn.GetP();

				vector<float> vRadii(iCount);
				m_pTree->neighbor_copy_radii_sq(&vRadii[0]);

				vector<long> vID(iCount*2);
				m_pTree->neighbor_copy_indices(&vID[0]);

				int i = 0 , j = 0;
				for(i=0;i<iCount;i+=2)
				{
					int iN1 = vID[i], iN2 = vID[i+1], iPos = 0;
					
					prob_t distSQ = vRadii[i/2];
					
					if(distSQ <= 0.0f) continue;
					
					//find sorted position in iN1 neighbor array to place element in
					for(j=0;j<vFound[iN1];j++)
					{
						if(distSQ < pnn[iN1][j].m_dist)
						{
							iPos = j;
							break;
						}
					}
					//do we need to shift elements after this one?
					if(iPos < iNNToFind - 1)
					{
						for(j=vFound[iN1]-1;j>iPos;j--)
							pnn[iN1][j]=pnn[iN1][j-1];
						pnn[iN1][iPos]=Neighbor(iN2,distSQ);
						if(vFound[iN1]+1<iNNToFind)vFound[iN1]++;
					}
					else if(iPos<iNNToFind)//no shift needed, just store
					{	
						pnn[iN1][iPos]=Neighbor(iN2,distSQ);
						if(vFound[iN1]+1<iNNToFind)vFound[iN1]++;
					}

					iPos = 0;
					for(j=0;j<vFound[iN2];j++)
					{
						if(distSQ < pnn[iN2][j].m_dist)
						{
							iPos = j;
							break;
						}
					}
					if(iPos < iNNToFind - 1)
					{
						for(j=vFound[iN2]-1;j>iPos;j--)
							pnn[iN2][j]=pnn[iN2][j-1];
						pnn[iN2][iPos]=Neighbor(iN1,distSQ);
						if(vFound[iN2]+1<iNNToFind)vFound[iN2]++;
					}
					else if(iPos<iNNToFind)	
					{
						pnn[iN2][iPos]=Neighbor(iN1,distSQ);
						if(vFound[iN2]+1<iNNToFind)vFound[iN2]++;
					}
				}
				return fRad;
			}
			//increase search radius
			fRad *= fRadFctr;
		}		
	}

	//get nearest neighbor as float*
	float* GetNearestNeighbor(float* p,bool bAllowZeroDist)
	{
		if(m_iNumElems == 1) return 0;
		
		const int iNumRads = 7;
		float pRads[7] = {3.0f,30.0f,150.0f,300.0f,600.0f,900.0f,1000.0f};

		int iIter = 0;
		float fRad = pRads[0];

		while(true)
		{
			m_pTree->search_center_radius_sq(p,fRad,1);
			int iCount = m_pTree->get_count();
			if(iCount)
			{
				vector<float> vRadii(iCount);
				m_pTree->copy_radii_sq(&vRadii[0]);

				int id = MinIdx(vRadii,bAllowZeroDist);
				if(id != -1)
				{
					vector<long> vID(iCount);
					m_pTree->copy_indices(&vID[0]);
					return &m_vData[vID[id]*m_iDims];
				}
			}
			//increase search radius
			if(iIter+1 >= iNumRads)
				fRad *= two;
			else
				fRad = pRads[++iIter];
		}		
	}

	float GetNearestRadiusSQ(float* p,vector<int>& vMap,int iID)
	{
		if(m_iNumElems == 1) return 0.0;
		
		extern prob_t gstartrad;
		float fRad = gstartrad;

		while(true)
		{
			m_pTree->search_center_radius_sq(p,fRad,1);
			int iCount = m_pTree->get_count();
			if(iCount)
			{
				vector<float> vRadii(iCount);
				m_pTree->copy_radii_sq(&vRadii[0]);

				float fm = FLT_MAX; 

				vector<long> vIndices(iCount);
				m_pTree->copy_indices(&vIndices[0]);

				bool bFound = false;

				int i = 0;
				for(i=0;i<iCount;i++)
				{
					if(vMap[vIndices[i]]==iID && vRadii[i]<= fm && vRadii[i]>0.0)
					{	bFound = true;
						fm = vRadii[i];
					}
				}
				if(bFound)
					return fm;
			}
			//increase search radius
			fRad *= two;
		}		
	}

	bool GetNearestNeighbor(float* p,bool bAllowZeroDist,Neighbor& n)
	{
		if(m_iNumElems == 1) return false;
		m_pTree->search_nn(p,bAllowZeroDist);
		float fm = 0.0;
		m_pTree->copy_radii_sq(&fm);
		long id;
		m_pTree->copy_indices(&id);
		n.m_dist = fm;
		n.m_id = id;
		return true;
	}

	float GetNearestRadiusSQ(float* p,bool bAllowZeroDist,bool bTest)
	{
		if(m_iNumElems == 1) return 0.0;		

		extern prob_t gstartrad;
		float fRad = gstartrad;

		if(bTest)
		{
#ifdef DO_TIMING
			extern MTimer oMT;	TimerInc oT(oMT);
#endif
			m_pTree->search_nn(p,bAllowZeroDist);
			float fm = 0.0;
			m_pTree->copy_radii_sq(&fm);
			return fm;
		}
		else
		{
#ifdef DO_TIMING
			extern MTimer oMF; TimerInc oT(oMF);
#endif
			while(true)
			{
				m_pTree->search_center_radius_sq(p,fRad,1);
				int iCount = m_pTree->get_count();
				if(iCount)
				{
					vector<float> vRadii(iCount);
					m_pTree->copy_radii_sq(&vRadii[0]);

					float fm = 0.0;

					if(bAllowZeroDist)
						return *std::min_element(vRadii.begin(),vRadii.end());
					else 
						fm=GZeroMinElem(vRadii);

					if(fm>0.0)
						return fm;
				}
				//increase search radius
				fRad *= two;
			}
		}
	}

	float GetNearestRadiusSQ(int i,bool bTest)
	{
		return GetNearestRadiusSQ(&m_vData[i*m_iDims],false,bTest);
	}

	//returns probability based on distance
	//of an arbitrary element in THIS distribution
	//to it's nearest neighbor in THIS distribution
	prob_t RProb(prob_t dRad)
	{
		return m_dTop / (m_dPiPow*dRad*m_dGamma);
	}

	//returns probability based on distance
	//of an arbitrary element in a DIFFERENT distribution
	//to it's nearest neighbor in THIS distribution
	prob_t RProbOther(prob_t dRad)
	{
		return (one/m_iNumElems) / (m_dPiPow*dRad*m_dGamma);
	}

	//returns probability of element i
	prob_t IProb(int i)
	{
		if(!m_pTree || i<0 || i>=m_iNumElems) return 0.0;

		if(1==m_iNumElems)return 1.0;

		return VProb(&m_vData[i*m_iDims]);
	}

	//returns probability of vector p
	//vector p must be in THIS distribution
	prob_t VProb(float* p)
	{
		if(!p || !m_pTree) return 0.0;

		if(1==m_iNumElems)return 1.0;

		prob_t dRad = sqrt(GetNearestRadiusSQ(p,false,false));

		return RProb(dRad);
	}

	//returns probability of vector p
	//vector p must be in DIFFERENT distribution
	prob_t VProbOther(float* p)
	{
		if(!p || !m_pTree) return 0.0;

		if(1==m_iNumElems)return 1.0;

		prob_t dRad = sqrt(GetNearestRadiusSQ(p,true,false));

		if(dRad == 0.0) return 0.0;

		return RProbOther(dRad);
	}

	char m_strMsg[1024];

	//entropy of distribution
	prob_t Entropy()
	{	//sprintf(m_strMsg,"Entropy sz=%d",m_iNumElems);
		//ScopedTimer S(m_strMsg);
		if(m_iNumElems<2) return 0.0;
		
		prob_t dEntrop = 0.0;
		prob_t dPiPowGamma = m_dPiPow*m_dGamma;

		int isz = m_iNumElems , i=0, iOffset = 0;
		
		for(i=0;i<isz;i++)
		{
			prob_t dDist = GetNearestRadiusSQ(&m_vData[iOffset],false,false);
			if(dDist<=0.0)continue;
			dDist = sqrt(dDist); 
			prob_t dProb = m_dTop / (dDist*dPiPowGamma);
			if(dProb<=0.0)continue;
			dEntrop += dProb * log2(dProb);
			iOffset += m_iDims;
		}
		return -dEntrop;
	}

	//vIDs specifies which cluster each element
	//belongs to. iClust specifies which cluster
	//to get entropy for
	prob_t Entropy(vector<int>& vIDs,int iClust)
	{
		prob_t dEntrop = 0.0;
		int isz = m_iNumElems , i=0;
		for(i=0;i<isz;i++)
		{
			if(vIDs[i]==iClust)
			{
				prob_t dProb = IProb(i);
				if(dProb==0.0) continue;
				dEntrop += dProb * log2(dProb);
			}
		}
		return -dEntrop;
	}

	//vIDs specifies which cluster each element
	//belongs to. iClust specifies which cluster
	//to get entropy for
	prob_t Entropy(vector<int>& vIDs,int iClust,int iNumElems)
	{	sprintf(m_strMsg,"Entropy c%d sz=%d totsz=%d",iClust,iNumElems,m_iNumElems);
	  //		ScopedTimer S(m_strMsg);
		if(iNumElems<2) return 0.0;

		const prob_t PI=3.14159265358979323846;
		prob_t dPiPow = pow((prob_t)PI,(prob_t)m_iDims/two);
		prob_t dTop = (one/(iNumElems-one));
		prob_t dGamma = Gamma(m_iDims/two+one);
		prob_t dPiPowGamma = dPiPow*dGamma;

		prob_t dEntrop = 0.0;
		int isz = m_iNumElems , i=0;
		for(i=0;i<isz;i++)
		{
			if(vIDs[i]==iClust)
			{
				prob_t dDist = GetNearestRadiusSQ(&m_vData[i*m_iDims],vIDs,iClust);
				if(dDist<=0.0)continue;
				dDist=sqrt(dDist);
				prob_t dProb = dTop / (dDist*dPiPowGamma);
				if(dProb<=0.0)continue;
				dEntrop += dProb * log2(dProb);
			}
		}
		return -dEntrop;
	}
};

//this is the continuous multidimensional probability version
inline void FillDistribs(vector<float>& vFloat,vector<KDTreeHist>& vDistribs,vector<KDTreeHist>& vCompDistribs,int iDistribs,vector<int>& vClustIDs,vector<int>& vCounts,int iDims,bool bGetComplements)
{
	vDistribs = vector< KDTreeHist >(iDistribs+1);

	if(bGetComplements)
		vCompDistribs = vector< KDTreeHist >(iDistribs+1);

	int iV = 0 , iC = 0;

	//full distribution
	vDistribs[iDistribs].SetData(iDims,&vFloat[0],vFloat.size()/iDims);

	int iTotalVs = vFloat.size() / iDims;

	for(iC=1;iC<iDistribs;iC++)
	{
		vector<float> vClustData(vCounts[iC]*iDims), vCompData;
		int iCompSize = iTotalVs - vCounts[iC];
		if(bGetComplements) vCompData = vector<float>(iCompSize*iDims);
		int j = 0 , k = 0;
		for(iV=0;iV<vClustIDs.size();iV++)
		{
			if(vClustIDs[iV] == iC)
			{
				int iD = 0;
				for(iD=0;iD<iDims;iD++)
				{
					vClustData[j++]=vFloat[iV*iDims+iD];
				}
			}
			else if(bGetComplements)
			{
				int iD = 0;
				for(iD=0;iD<iDims;iD++)
				{
					vCompData[k++]=vFloat[iV*iDims+iD];
				}
			}
		}
		vDistribs[iC].SetData(iDims,&vClustData[0],vCounts[iC]);
		if(bGetComplements)
			vCompDistribs[iC].SetData(iDims,&vCompData[0],iCompSize);
	}
}

//this is the continuous multidimensional probability version
void FillDistribs(vector<float>& vFloat,vector<KDTreeHist>& vDistribs,vector<KDTreeHist>& vCompDistribs,int iDistribs,vector<int>& vClustIDs,vector<int>& vCounts,int iDims,A2D<int>& vBestDims,int iBestDims);
//get full 'background' distribution , containing all spikes but using dimensions specified in pBestDims
void GetFullBGDistrib(vector<float>& vFloat,KDTreeHist& oTree,int iDims,int* pBestDims,int iBestDims);

// bool RandAssign(CVerxStack& DataStack,CCluster& MainClusters,int iClusts,int which);

#endif
