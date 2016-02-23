//////////////////////////////////////////////////////////////////////
// Cluster.cpp: implementation of the CCluster class.
//
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

//////////////////////////////////////////////////////////////////////

#include <float.h>
#include "Vertex.h"
#include "Cluster.h"
#include "IsolationDist.h"
#include <algorithm>
#include <map>

#include "Log.h"
#include "hr_time.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// CCluster
CCluster::CCluster() { 
  m_iNumClusts = 0;	
};

//////////////////////////////////////////////////////////////////////
void CCluster::Clear() {
  m_iNumClusts = 0;	
}

int CCluster::GetNumClusts() {  
  return m_iNumClusts;
}

bool ClustIsolationDist(CCluster& Clusters,CVerxStack& oS,vector<double>& vFloat,vector<int>& vClustIDs,vector<int>& vClustCounts,int iClusts,int iRows,int iCols,bool bTime)
{	try{
	char msg[1024];
	sprintf(msg,"Computing isolation distance for clusters");

	Write2Log("Computing isolation distance for clusters");
	int iC = 0;
	if(Clusters.m_vInfo.size() < iClusts+1)	Clusters.m_vInfo.resize(iClusts+1);
	CStopWatch oTimer;	
	vector<float> vRange(iCols),vMin(iCols);
	for(iC=1;iC<=iClusts;iC++)
	{	sprintf(msg,"Computing isolation distance for cluster %d",iC); Write2Log(msg);	        
		if(vClustCounts[iC]<2)
		{	Write2Log("Clust too small for isolation distance : %d %d , skipping",iC,vClustCounts[iC]);
			continue;
		}
		float fDist = 0.0;
		try
		{	
		  if(bTime) oTimer.startTimer();
		  fDist=IsolationDist(vClustIDs,iC,vFloat,iRows,iCols);
		  if(bTime) oTimer.stopTimer();
		}
		catch(...)
		{	Write2Log("Exception computing isolation dist on clust=%d sz=%d",iC,vClustCounts[iC]);
		}
		if(bTime) Clusters.m_vInfo[iC].m_dIsoDTime = oTimer.getElapsedTime();
		Clusters.m_vInfo[iC].m_dIsolationDist = fDist;
		Write2Log("Isolation dist of clust %d is %.6f",iC,fDist);
	}
	return true;
	}catch(...){
		Write2Log("Outer exception in ClustIsolationDist!!");
		return false;
	}
}

bool ClustLRatio(CCluster& Clusters,CVerxStack& oS,vector<double>& vFloat,vector<int>& vClustIDs,vector<int>& vClustCounts,int iClusts,int iRows,int iCols,bool bTime)
{	try{
	char msg[1024];
	sprintf(msg,"Computing L-Ratio for clusters");
	Write2Log("Computing L-Ratio for clusters");
	int iC = 0;
	if(Clusters.m_vInfo.size() < iClusts+1)	Clusters.m_vInfo.resize(iClusts+1);
	CStopWatch oTimer;
	vector<float> vRange(iCols),vMin(iCols);
	for(iC=1;iC<=iClusts;iC++)
	{	sprintf(msg,"Computing L-Ratio for cluster %d",iC);	
		Write2Log(msg);
		float fLRatio = -1.0;
		try
		{
		  if(bTime) oTimer.startTimer();
		  fLRatio=LRatio(vClustIDs,iC,vFloat,iRows,iCols);
		  if(bTime) oTimer.stopTimer();
		}
		catch(...)
		{ Write2Log("Exception computing L-Ratio on clust=%d sz=%d",iC,vClustCounts[iC]);
		}
		if(bTime) Clusters.m_vInfo[iC].m_dLRatTime = oTimer.getElapsedTime();
		Clusters.m_vInfo[iC].m_dLRatio = fLRatio;
		Write2Log("L-Ratio of clust %d is %.6f",iC,fLRatio);
	}
	return true;
	}catch(...){
		Write2Log("Outer exception in ClustLRatio!!");
		return false;
	}
}
