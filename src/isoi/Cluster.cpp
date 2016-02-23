// $Id: Cluster.cpp,v 1.6 2011/09/02 16:12:12 samn Exp $ 
// Cluster.cpp: implementation of the CCluster class.
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

#include <float.h>
#include "Vertex.h"
#include "Cluster.h"
#include "Hist.h"
#include <algorithm>
#include <map>
#include "KDTree.h"
#include "Log.h"
#include "InfoT.h"
#include "hr_time.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// CCluster
CCluster::CCluster() { 
  m_iNumClusts = 0;	
  m_oCQO.Default();
};

//////////////////////////////////////////////////////////////////////
void CCluster::Clear() {
  m_iNumClusts = 0;	
}

int CCluster::GetNumClusts() {  
  return m_iNumClusts;
}

prob_t gstartrad = 0.1;

bool CalcClustInfo(CVerxStack& DataStack,bool bSymmetric,A2D<Neighbor>& vnn,vector<float>& vFloat,
		   KDTreeHist& oTree,vector<float>& vRange,bool bWC,bool bTime) {
  extern int iSQRTCalls; iSQRTCalls = 0;
  CCluster& Clusters = DataStack.m_oClusters;
  CStopWatch oTimer;
  try {
    char msg[1024]={0};

    DataStack.CalcMinMax();
		
    int iClusts = Clusters.GetNumClusts() , iC=1;
    int iDims=    DataStack.m_iNumDim, iD=0;

    vector<int> vCounts(iClusts+1);
    DataStack.GetCounts(vCounts,iClusts);
		
    //for info estimation use Peak-V, even though it is not used in auto-clustering
    vector<int> vZeroes(iDims);
		
    if(Clusters.m_vInfo.size()<iClusts+1)
      Clusters.m_vInfo = vector< ClusterInfo >(iClusts+1);
    for(iC=1;iC<=iClusts;iC++) {
      Clusters.m_vInfo[iC].m_iMyID = iC;
      Clusters.m_vInfo[iC].m_iSz = vCounts[iC];
    }

    //multidimensional distributions - continuous
    //resistor average of all points in each cluster against everything OTHER than that cluster points
    int iRows=0,iCols=0;
    DataStack.GetFloatV(iRows,iCols,vFloat,vRange);

    vector<int> vClustIDs;
    DataStack.GetClustIDs(vClustIDs);

    vector< KDTreeHist > vDistribs, vCompDistribs;

    const bool bFindBestDims = Clusters.m_oCQO.m_bFindBestDims;
    A2D<int>& vBestDims = Clusters.m_vBestDims;
    A2D<prob_t>& v1DKLDivs = Clusters.m_v1DKLDivs;
    A2D<KLD2D>& v2DKLDivs = Clusters.m_v2DKLDivs;
    int iBestDims = Clusters.m_oCQO.m_iBestDims;
    const bool b1D = Clusters.m_oCQO.m_b1DBest;
    const bool bUseDefDims = false; // Clusters.m_oCQO.m_bUseDefDims;
    if(bFindBestDims) {
      sprintf(msg,"Finding best %d dimensions for each of the %d clusters.",iBestDims,iClusts);
      Write2Log(msg);
      Clusters.m_oCQO.m_bFastKLD = false;
      if(b1D){	//find best dims using 1D slices
		FindBest1DDims(vFloat,iClusts,iCols,iBestDims,vCounts,vClustIDs,vBestDims,v1DKLDivs);
		//printout for debuggering
		Write2Log("v1DKLDivs:\n");	WriteMat2Log(v1DKLDivs);
		Write2Log("\nvBestDims");	WriteMat2Log(vBestDims);
      }
      else	//find best dims using 2D slices
	        FindBest2DDims(vFloat,vRange,iClusts,iCols,iBestDims,vCounts,vClustIDs,vBestDims,v2DKLDivs,bWC,bTime,Clusters.m_vInfo);
      FillDistribs(vFloat,vDistribs,vCompDistribs,iClusts+1,vClustIDs,vCounts,iCols,vBestDims,iBestDims);
    }
    else if(bUseDefDims) { //use 8 dimensions of T1-4Peak,T1-4Energy
      vBestDims.Init(iClusts+1,8);
      iBestDims = 8;
      for(iC=1;iC<=iClusts;iC++) {
		for(iD=0;iD<=3;iD++) vBestDims[iC][iD]=iD; //T1-4Peak
		for(iD=16;iD<=19;iD++) vBestDims[iC][iD-12]=iD; //T1-4Energy
      }
      int x,y; for(y=1;y<vBestDims.Rows() && y<2;y++) {
	  	Write2Log("\ndefault dims:\t",y);
	  	for(x=0;x<vBestDims.Cols();x++) Write2Log("index=%d\t",vBestDims[y][x]);
	  	Write2Log("\n\n");
      }
      FillDistribs(vFloat,vDistribs,vCompDistribs,iClusts+1,vClustIDs,vCounts,iCols,vBestDims,iBestDims);
      Write2Log("filled distribs");
    }
    else //use all available dimensions
      FillDistribs(vFloat,vDistribs,vCompDistribs,iClusts+1,vClustIDs,vCounts,iCols,true);

    bool bFast = Clusters.m_oCQO.m_bFastKLD && (Clusters.m_oCQO.m_bDoKLDiv || Clusters.m_oCQO.m_bInterKLD);
    int iNNToFind = 1; 
		
    if(Clusters.m_oCQO.m_bInterKLD || Clusters.m_oCQO.m_bDoKLDiv)
      Write2Log("Starting CalcClustInfo bFast=%s NumElems=%d NumDims=%d",bFast?"true":"false",iRows,iCols);
    Neighbor** pnn = 0; //neighbors , only used in fast mode
    vector<int> vFound(iRows,1);
    if(bFast) {
      if(!vnn.Init(iRows,iNNToFind)) {
		Write2Log("CIThread: Out of memory, couldn't alloc %d by %d Neighbors",iRows,iNNToFind);
		return false;
      }
      pnn = vnn.GetP();
      //find 1st nearest neighbors for each vertex
      oTree.SetData(iCols,&vFloat[0],iRows);
      sprintf(msg,"Finding initial %d nearest neighbors...",iNNToFind);			
      int k;
      for(k=0;k<iRows;k++) {			  	
	if(k%1000==0)
	  fprintf(stderr,"Finding nearest neighbor for vector %d",k);
	oTree.GetNearestNeighbor(oTree[k],false,pnn[k][0]);
      }
    }

    if(bFast) {
    	NeighborClustCount(pnn,vClustIDs,Clusters.m_vInfo,iNNToFind,iClusts);
      	int s = 0;
      	for(s=1;s<=iClusts;s++) {		
			Write2Log("clust %d avg prct %d-nearest-neighbors in same clust = %.4f",s,iNNToFind,Clusters.m_vInfo[s].m_fPrctKNNInClust);
      	}
    }		
    vector<prob_t> vcInf(iClusts+1);//kldiv resistor avg. cluster vs background
    CStopWatch oTimer;
    if(Clusters.m_oCQO.m_bDoKLDiv) for(iC=1;iC<=iClusts;iC++) {
	try {			  
	  fprintf(stderr,"Calculating IsoI_BG for cluster %d of %d\n",iC,iClusts);
	} catch(...) {	
	  Write2Log("Exception in progress control Cluster=%d Count=%d",iC,vCounts[iC]);
	}
	prob_t kldiv=0.0f;	      
	if(vCounts[iC]) {			
	  try {
	    if(bTime) oTimer.startTimer();
	    if(bFast) {	
	    	kldiv = FastKLDivSymPNOTP(vDistribs[iC],vCompDistribs[iC],pnn,iC,vClustIDs,iNNToFind,vFound);
	    }
	    else {					
	    	kldiv = bSymmetric ? 
			KLDivSym(vDistribs[iC],vCompDistribs[iC]) : KLDiv(vDistribs[iC],vCompDistribs[iC]);
	    }
	    if(bTime) { // save time to get IsoI_BG
	      oTimer.stopTimer();
	      Clusters.m_vInfo[iC].m_dBGTime = oTimer.getElapsedTime();
	    }
	    Write2Log(bSymmetric?"sym. kldiv=%.4f":"kldiv=%.4f",kldiv);
	  } catch(...) {
	  	Write2Log("Exception in KLDivSym Cluster=%d Count=%d",iC,vCounts[iC]);
	  }
	  if(kldiv<0.0) Write2Log("kldiv<0.0 Cluster=%d Count=%d",iC,vCounts[iC]);
	}
	else
	  Write2Log("Empty Cluster!!! Cluster=%d Count=0",iC);
	try {
	  	//resistor average of kldiv from cluster iC to background
	  	Clusters.m_vInfo[iC].m_fBGInfoGain = kldiv;
	}
	catch(...) {
	  	Write2Log("exception in assigning Clusters.m_vInfo Cluster=%d Count=%d",iC,vCounts[iC]);
	}
   }
    if(Clusters.m_oCQO.m_bInterKLD) //kldiv cluster to other clusters
      InterClustKLD(Clusters,vDistribs,vClustIDs,vCounts,iClusts,bFast,pnn,0,iNNToFind,vFound,vFloat,iRows,iCols,bTime);

    if(Clusters.m_oCQO.m_bDoKLDiv) for(iC=1;iC<=iClusts;iC++)
				     Write2Log("clust=%d count=%d kldiv_total=%.4f kldiv_bg=%.4f kldiv_inter=%.4f closest_clust=%d",
					       iC,vCounts[iC],
					       Clusters.m_vInfo[iC].m_fBGInfoGain+Clusters.m_vInfo[iC].m_fInterClustGain,
					       Clusters.m_vInfo[iC].m_fBGInfoGain,
					       Clusters.m_vInfo[iC].m_fInterClustGain,
					       Clusters.m_vInfo[iC].m_iClosestID);

    if((Clusters.m_oCQO.m_bDoKLDiv || Clusters.m_oCQO.m_bInterKLD))
      Write2Log("Finished CalcClustInfo bFast=%d bDoKLDiv=%d bInter=%d NumElems=%d NumDims=%d iSQRTCalls=%d bFindBestDims=%d iBestDims=%d",
		bFast?1:0,Clusters.m_oCQO.m_bDoKLDiv?1:0,Clusters.m_oCQO.m_bInterKLD?1:0,iRows,iCols,iSQRTCalls,Clusters.m_oCQO.m_bFindBestDims?1:0,Clusters.m_oCQO.m_iBestDims);
		
    return true;
  } 
  catch(...) {
    Write2Log("CIThread: Outermost Exception!");
    return false;
  }
}
