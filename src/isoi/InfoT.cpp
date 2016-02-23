// $Id: InfoT.cpp,v 1.7 2011/09/02 16:10:57 samn Exp $ 
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

#include "InfoT.h"
#include "WCMath.h"
#include "combination.h"
#include "hr_time.h"

//#define LOG_TEST 1

#ifdef DO_TIMING
 MTimer oMT("getnnrsq btest=true");
 MTimer oMF("getnnrsq btest=false");
#endif

#define SLOWER_KD 1
#define FASTER_KD 2
#define BOTH_KD 4
#define KD_MODE FASTER_KD

bool NeighborClustCount(Neighbor** vNeighbors,vector<int>& vClustIDs,vector<ClusterInfo>& vPrct,int iNNToFind,int iClusts)
{
	int iSz = vClustIDs.size();
	if(vPrct.size()!=iClusts+1)
		return false;
	int i;
	vector<int> vCounts(iClusts+1);
	for(i=1;i<=iClusts;i++)
		vCounts[i]=count(vClustIDs.begin(),vClustIDs.end(),i);
	vector< vector<prob_t> > vprct(iClusts+1);
	for(i=1;i<=iClusts;i++)
		vprct[i]=vector<prob_t>(vCounts[i],0.0);
	vector<int> vIndex(iClusts+1,0);
	for(i=0;i<iSz;i++)
	{
		int j, iClust = vClustIDs[i];
		if(!iClust)continue;//skip background vectors
		int idx = vIndex[iClust]++;
		for(j=0;j<iNNToFind;j++)
		{
			if(vClustIDs[vNeighbors[i][j].m_id]==iClust)
				vprct[iClust][idx]+=1.0f;
		}
		vprct[iClust][idx] /= iNNToFind;
	}
	for(i=1;i<=iClusts;i++)
		vPrct[i].m_fPrctKNNInClust=Avg(vprct[i]);
	return true;
}

prob_t FastKLDivPQ(KDTreeHist& p,KDTreeHist& q,Neighbor** vNeighbors,int iClustIDP,int iClustIDQ,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount)
{
	int iFastP = 0, iFastQ = 0, iSlowP = 0 , iSlowQ = 0;

#ifdef LOG_TEST
	double kldiv = 1.0;
	double kldivQ=1.0,kldivP=1.0;
#else
	prob_t kldiv = 0.0;
#endif
	int isz = vClustIDs.size() , iV = 0, iT = 0, iNN = 0, iLessP = 0, iLessQ = 0;

	for(iV=0;iV<isz;iV++)
	{
		if(vClustIDs[iV]!=iClustIDP)	//skip vectors not in cluster p
			continue;

		prob_t dDistP = 0.0 , dDistQ = 0.0; // squared distances
		
		Neighbor* vnn = vNeighbors[iV]; //these are the neighbors of an element in cluster p
		int nsz = vNCount[iV]; //iNNToFind;//vnn.size();
		bool bFoundP = false, bFoundQ = false;
		for(iNN=0;iNN<nsz;iNN++)
		{	Neighbor& oN = vnn[iNN];
			if(!bFoundP && vClustIDs[oN.m_id]==iClustIDP)
			{	if(oN.m_dist>0.0)
				{	//found a different neighbor in the same cluster p
					dDistP = oN.m_dist;
					iFastP++;
					bFoundP = true;
				}
			}
			else if(!bFoundQ && vClustIDs[oN.m_id]==iClustIDQ)
			{	//found a neighbor in the other cluster q
				dDistQ = oN.m_dist;
				iFastQ++;
				bFoundQ = true;
			}
			if(bFoundP && bFoundQ) break;	//found neighbors so break
		}
		if(!bFoundP)
		{	//use slow method of searching KD-tree
#if KD_MODE == BOTH_KD // to make sure get same results
			dDistP = p.GetNearestRadiusSQ(iT,false);
			float dDistP2 = p.GetNearestRadiusSQ(iT,true);
			if(dDistP!=dDistP2)
			{
				Write2Log("not eq %.6f %.6f",dDistP,dDistP2);
			}
#elif KD_MODE == SLOWER_KD // in case have to revert for some unkown reason...
			dDistP = p.GetNearestRadiusSQ(iT,false);
#else
			dDistP = p.GetNearestRadiusSQ(iT,true);
#endif
			iSlowP++;
		}
		if(!bFoundQ)
		{	//use slow method of searching KD-tree
#if KD_MODE == BOTH_KD // to make sure get same results
			dDistQ = q.GetNearestRadiusSQ(p[iT],true,false);//true==alow zero distance,since diff cluster,can have same exact point with 0 dist
			float dDistQ2 = q.GetNearestRadiusSQ(p[iT],true,true);
			if(dDistQ2!=dDistQ)
			{
				Write2Log("not eq %.6f %.6f",dDistQ,dDistQ2);
			}
#elif KD_MODE == SLOWER_KD
			dDistQ = q.GetNearestRadiusSQ(p[iT],true,false);//true==alow zero distance,since diff cluster,can have same exact point with 0 dist
#else
			dDistQ = q.GetNearestRadiusSQ(p[iT],true,true);//true==alow zero distance,since diff cluster,can have same exact point with 0 dist
#endif
			iSlowQ++;
		}

		if(dDistP>0.0 && dDistQ>0.0)
		{	//avoid exceptions of log(0)
			//if(dDistP<dDistQ)iLessP++; else iLessQ++;
#ifdef LOG_TEST
			//kldiv *= dDistQ / dDistP;
	//		kldivQ *= dDistQ;
	//		kldivP *= dDistP;
			kldiv += log( dDistQ / dDistP );
#else
			kldiv += log2( dDistQ / dDistP ) / 2.0;
#endif
		}
		
		iT++;	//increment index into cluster p's KD-tree
	}
	//finish the calculation
#ifdef LOG_TEST
	//kldiv = log2( kldiv ) / 2.0;
	//kldiv = log2( kldivQ / kldivP ) / 2.0;
	kldiv /= ( log(2.0) * 2.0 );
#endif
	kldiv *= p.NumDims() / ((prob_t) p.NumElems() );
	kldiv += log2( (prob_t)q.NumElems() / (p.NumElems()-1.0 ) );
	//write some stats
	Write2Log("FastKLDivPQ:kldiv=%.4f iSlow=%d iSlowP=%d iSlowQ=%d iFast=%d iFastP=%d iFastQ=%d",kldiv,iSlowP+iSlowQ,iSlowP,iSlowQ,iFastP+iFastQ,iFastP,iFastQ);
	
	return kldiv;
}

prob_t FastKLDivPNOTP(KDTreeHist& p,KDTreeHist& notp,Neighbor** vNeighbors,int iClustIDP,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount)
{
	int iFastP = 0, iFastNotP = 0, iSlowP = 0 , iSlowNotP = 0;

#ifdef LOG_TEST
	double kldiv = 1.0;
	double kldivP=1.0,kldivNOTP=1.0;
#else
	prob_t kldiv = 0.0;
#endif
	int isz = vClustIDs.size() , iV = 0, iT = 0, iNN = 0, iLessP = 0;

	for(iV=0;iV<isz;iV++)
	{
		if(vClustIDs[iV]!=iClustIDP)	//skip vectors not in cluster p
			continue;

		prob_t dDistP = 0.0 , dDistNotP = 0.0; // squared distances
		
		Neighbor* vnn = vNeighbors[iV]; //these are the neighbors of an element in cluster p
		int nsz = vNCount[iV]; //iNNToFind;//vnn.size();
		bool bFoundP = false, bFoundNotP = false;
		for(iNN=0;iNN<nsz;iNN++)
		{	Neighbor& oN = vnn[iNN];
			if(!bFoundP && vClustIDs[oN.m_id]==iClustIDP)
			{	if(oN.m_dist>0.0)
				{	//found a different neighbor in the same cluster p
					dDistP = oN.m_dist;
					iFastP++;
					bFoundP = true;
				}
			}
			else if(!bFoundNotP && vClustIDs[oN.m_id]!=iClustIDP)
			{	//found a neighbor in the other cluster q
				dDistNotP = oN.m_dist;
				iFastNotP++;
				bFoundNotP = true;
			}
			if(bFoundP && bFoundNotP) break;	//found neighbors so break
		}
		if(!bFoundP)
		{	//use slow method of searching KD-tree
#if KD_MODE == BOTH_KD // to make sure get same results
			dDistP = p.GetNearestRadiusSQ(iT,false);
			float dDistP2 = p.GetNearestRadiusSQ(iT,true);
			if(dDistP!=dDistP2)
			{
				Write2Log("not eq %.6f %.6f",dDistP,dDistP2);
			}
#elif KD_MODE == SLOWER_KD
			dDistP = p.GetNearestRadiusSQ(iT,false);
#else
			dDistP = p.GetNearestRadiusSQ(iT,true);
#endif
			iSlowP++;
		}
		if(!bFoundNotP)
		{	//use slow method of searching KD-tree
#if KD_MODE == BOTH_KD // to make sure get same results
			dDistNotP = notp.GetNearestRadiusSQ(p[iT],true,false);//true==alow zero distance,since diff cluster,can have same exact point with 0 dist
			float dDistNotP2 = notp.GetNearestRadiusSQ(p[iT],true,true);
			if(dDistNotP!=dDistNotP2)
			{
				Write2Log("not eq %.6f %.6f",dDistNotP,dDistNotP2);
			}
#elif KD_MODE == SLOWER_KD
			dDistNotP = notp.GetNearestRadiusSQ(p[iT],true,false);//true==alow zero distance,since diff cluster,can have same exact point with 0 dist
#else
			dDistNotP = notp.GetNearestRadiusSQ(p[iT],true,true);//true==alow zero distance,since diff cluster,can have same exact point with 0 dist
#endif
			iSlowNotP++;
		}

		if(dDistP>0.0 && dDistNotP>0.0)
		{	//avoid exceptions of log(0)
			//if(dDistP<dDistNotP)iLessP++; else iLessNotP++;
#ifdef LOG_TEST
			//kldiv *= dDistNotP / dDistP;
			//kldivP *= dDistP;
			//kldivNOTP *= dDistNotP;
			kldiv += log ( dDistNotP / dDistP );
#else
			kldiv += log2( dDistNotP / dDistP ) / 2.0;
#endif
		}
		
		iT++;	//increment index into cluster p's KD-tree
	}
	//finish the calculation
#ifdef LOG_TEST
	//kldiv = log2( kldiv ) / 2.0;
	//kldiv = log2( kldivNOTP / kldivP ) / 2.0;
	kldiv /= ( log(2.0) * 2.0 );
#endif
	kldiv *= p.NumDims() / ((prob_t) p.NumElems() );
	kldiv += log2( (prob_t)notp.NumElems() / (p.NumElems()-1.0 ) );
	//write some stats
	Write2Log("FastKLDivPNOTP:kldiv=%.4f iSlow=%d iSlowP=%d iSlowNotP=%d iFast=%d iFastP=%d iFastNotP=%d",kldiv,iSlowP+iSlowNotP,iSlowP,iSlowNotP,iFastP+iFastNotP,iFastP,iFastNotP);
	
	return kldiv;
}

//notp is all points not in p
//this version gets the distance from notp to p
prob_t FastKLDivNOTPP(KDTreeHist& p,KDTreeHist& notp,Neighbor** vNeighbors,int iClustIDP,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount)
{
	int iFastP = 0, iFastNP = 0, iSlowP = 0, iSlowNP = 0;

#ifdef LOG_TEST
	double kldiv = 1.0;
	double kldivNOTP=1.0,kldivP=1.0;
#else
	prob_t kldiv = 0.0;
#endif
	int isz = vClustIDs.size() , iV = 0, iT = 0, iNN = 0, iLessP = 0;

	for(iV=0;iV<isz;iV++)
	{
		if(vClustIDs[iV]==iClustIDP) //skip all points in p
			continue;

		prob_t dDistP = 0.0 , dDistNotP = 0.0;
		
		Neighbor* vnn = vNeighbors[iV];
		int nsz = vNCount[iV]; //iNNToFind;
		bool bFoundP = false, bFoundNotP = false;
		for(iNN=0;iNN<nsz;iNN++)
		{	Neighbor& oN = vnn[iNN];
			if(!bFoundP && vClustIDs[oN.m_id]==iClustIDP)
			{
				dDistP = oN.m_dist;
				iFastP++;
				bFoundP = true;
			}
			else if(!bFoundNotP && vClustIDs[oN.m_id]!=iClustIDP)
			{
				if(oN.m_dist>0.0)
				{
					dDistNotP = oN.m_dist;
					iFastNP++;
					bFoundNotP = true;
				}
			}
			if(bFoundP && bFoundNotP) break;
		}
		if(!bFoundP)
		{
#if KD_MODE == BOTH_KD // to make sure get same results
			dDistP = p.GetNearestRadiusSQ(notp[iT],true,false);
			float dDistP2 = p.GetNearestRadiusSQ(notp[iT],true,true);
			if(dDistP!=dDistP2)
			{
				Write2Log("not eq %.6f %.6f",dDistP,dDistP2);
			}
#elif KD_MODE == SLOWER_KD
			dDistP = p.GetNearestRadiusSQ(notp[iT],true,false);
#else
			dDistP = p.GetNearestRadiusSQ(notp[iT],true,true);
#endif
			iSlowP++;
		}
		if(!bFoundNotP)
		{
#if KD_MODE == BOTH_KD // to make sure get same results
			dDistNotP = notp.GetNearestRadiusSQ(iT,false);
			float dDistNotP2 = notp.GetNearestRadiusSQ(iT,true);
			if(dDistNotP!=dDistNotP2)
			{
				Write2Log("not eq %.6f %.6f",dDistNotP,dDistNotP2);
			}
#elif KD_MODE == SLOWER_KD
			dDistNotP = notp.GetNearestRadiusSQ(iT,false);
#else
			dDistNotP = notp.GetNearestRadiusSQ(iT,true);
#endif
			iSlowNP++;
		}

		if(dDistP>0.0 && dDistNotP>0.0)
		{
			//if(dDistNotP<dDistP)iLessNotP++; else iLessP++;
#ifdef LOG_TEST
			//kldiv *= dDistP / dDistNotP;
			//kldivP *= dDistP;
			//kldivNOTP *= dDistNotP;
			kldiv += log ( dDistP / dDistNotP );
#else
			kldiv += log2( dDistP / dDistNotP ) / 2.0;
#endif
		}
		
		iT++;
	}
#ifdef LOG_TEST
	//kldiv = log2( kldiv ) / 2.0;
	//kldiv = log2( kldivP / kldivNOTP ) / 2.0;
	kldiv /= ( log(2.0) * 2.0 );
#endif
	kldiv *= notp.NumDims() / ((prob_t) notp.NumElems() );
	kldiv += log2( (prob_t)p.NumElems() / (notp.NumElems()-1.0 ) );

	Write2Log("FastKLDivNOTPP:kldiv=%.4f iSlow=%d iSlowP=%d iSlowNP=%d iFast=%d iFastP=%d iFastNP=%d",kldiv,iSlowP+iSlowNP,iSlowP,iSlowNP,iFastP+iFastNP,iFastP,iFastNP);
	
	return kldiv;
}

//symmetric kldiv 
//(dpq * dqp ) / (dpq + dqp)
//dpq == distance from p to q
//dqp == distance from q to p
//both q and p are actual clusters with valid IDs!!
prob_t FastKLDivSymPQ(KDTreeHist& p,KDTreeHist& q,Neighbor** vNeighbors,int iClustIDP,int iClustIDQ,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount)
{
	prob_t dpq = FastKLDivPQ(p,q,vNeighbors,iClustIDP,iClustIDQ,vClustIDs,iNNToFind,vNCount), //dist from p to q
		   dqp = FastKLDivPQ(q,p,vNeighbors,iClustIDQ,iClustIDP,vClustIDs,iNNToFind,vNCount); //dist from q to p
	if(!dpq && !dqp) return 0.0;
	return dpq*dqp / (dpq+dqp);
}

//symmetric kldiv
//(dpnotp * dnotpp) / (dpnotp + dnotpp)
//dpnotp == distance from p to not p
//dnotpp == distance from not p to p
//p == actual cluster with valid ID
//notp == complement of p (all vectors not in p), without an actual ID 
prob_t FastKLDivSymPNOTP(KDTreeHist& p,KDTreeHist& notp,Neighbor** pNeighbors,int iClustIDP,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount)
{
	prob_t dpq = FastKLDivPNOTP(p,notp,pNeighbors,iClustIDP,vClustIDs,iNNToFind,vNCount),  //dist from p to notp
		   dqp = FastKLDivNOTPP(p,notp,pNeighbors,iClustIDP,vClustIDs,iNNToFind,vNCount); //dist from notp to p
	if(!dpq && !dqp) return 0.0;
	return dpq*dqp / (dpq+dqp);
}

//make tree using dimension indices stored in pBestDims
//iRows is # of elements in full vFloat matrix, iCols is full # of dimensions in vFloat
void FillTree(vector<float>& vFloat,int iRows,int iCols,int iClustID,int iClustSz,vector<int>& vClustIDs,int* pBestDims,int iBestDims,KDTreeHist& oT,vector<float>& vClustData)
{
	vClustData.resize(iClustSz*iBestDims);
	int iV = 0 , idx = 0;
	for(;iV<iRows;iV++)
	{	if(vClustIDs[iV]!=iClustID) continue;
		int iD = 0;
		for(;iD<iBestDims;iD++)
			vClustData[idx++]=vFloat[iV*iCols+pBestDims[iD]];
	}
	oT.SetData(iBestDims,&vClustData[0],iClustSz);
}

//computes kldiv between each pair of actual clusters and adds the min to CCluster::m_vInfo
bool InterClustKLD(CCluster& Clusters,vector<KDTreeHist>& vDistribs,vector<int>& vClustIDs,vector<int>& vClustCounts,int iClusts,bool bFast,Neighbor** vnn,int WhichDraw,int iNNToFind,vector<int>& vNCount,vector<float>& vFloat,int iRows,int iCols,bool bTime) {
  //kldiv cluster to other clusters
  int iC1 = 0, iC2 = 0, iTot = ((iClusts-1)*(iClusts-1)+iClusts-1)/2, iCurr = 0;
  CStopWatch oTimer;
  try {	
    char msg[2048]; 
    vector< vector<prob_t> > vcInfInter(iClusts+1, vector<prob_t>(iClusts+1));		
    Write2Log("Calculating inter-cluster KLDiv");
    if(Clusters.m_oCQO.m_bFindBestDims)	{ //compute inter-clust kldiv using best dimensions
      iTot = iClusts*iClusts;		//distances are not symmetrical since use different dimensions
      for(iC1=1;iC1<=iClusts;iC1++) {
	if(bTime) Clusters.m_vInfo[iC1].m_dNNTime = 0.0;
	for(iC2=1;iC2<=iClusts;iC2++,iCurr++) {	
	  if(iC2==iC1) continue;
	  sprintf(msg,"Calculating IsoI_NN(%d,%d)",iC1,iC2);
	  fprintf(stderr,"%s\n",msg);					
	  Write2Log(msg);
	  KDTreeHist oT; vector<float> vTmpData;	//make temporary tree
	  if(bTime) oTimer.startTimer();
	  FillTree(vFloat,iRows,iCols,iC2,vClustCounts[iC2],vClustIDs,Clusters.m_vBestDims[iC1],Clusters.m_oCQO.m_iBestDims,oT,vTmpData);
	  vcInfInter[iC1][iC2]=KLDivSym(vDistribs[iC1],oT);
	  if(bTime) { // save IsoI_NN time, for all pairs for this cluster
	    oTimer.stopTimer();
	    Clusters.m_vInfo[iC1].m_dNNTime += oTimer.getElapsedTime();	    
	  }
	}
      }
    } else {	
      for(iC1=1;iC1<=iClusts;iC1++) {
	if(bTime) Clusters.m_vInfo[iC1].m_dNNTime = 0.0;
	for(iC2=iC1+1;iC2<=iClusts;iC2++,iCurr++) {
	  sprintf(msg,"Calculating IsoI_NN(%d,%d)",iC1,iC2);
	  fprintf(stderr,"%s\n",msg);
	  if(bTime) oTimer.startTimer();
	  if(bFast)
	    vcInfInter[iC1][iC2]=vcInfInter[iC2][iC1]=FastKLDivSymPQ(vDistribs[iC1],vDistribs[iC2],vnn,iC1,iC2,vClustIDs,iNNToFind,vNCount);
	  else
	    vcInfInter[iC1][iC2]=vcInfInter[iC2][iC1]=KLDivSym(vDistribs[iC1],vDistribs[iC2]);
	  if(bTime) { // save IsoI_NN time, for all pairs for this cluster
	    oTimer.stopTimer();
	    Clusters.m_vInfo[iC1].m_dNNTime += oTimer.getElapsedTime();	    
	  }
	}
      }
    }
    for(iC1=1;iC1<=iClusts;iC1++) {	
      prob_t min_int = iClusts>1 ? INF : 0.0; //*INF;
      int min_ind = 0;
      if(iClusts>1) 
	for(iC2=1;iC2<=iClusts;iC2++) {
	  if(iC1==iC2 || vClustCounts[iC2]<2)continue;
	  prob_t tmpK = vcInfInter[iC1][iC2];
	  if(tmpK<min_int) {	
	    min_int=tmpK; 
	    min_ind=iC2; 
	  }
	}
      Clusters.m_vInfo[iC1].m_fInterClustGain = min_int;
      Clusters.m_vInfo[iC1].m_iClosestID = min_ind;
      Write2Log("Nearest kldiv from clust %d to %d is %.6f",iC1,min_ind,min_int);
    }
    string strTab("\ninter clust kldiv table\n");char strTmp[1024];//write inter-cluster kldiv table to log for inspection...
    for(iC1=1;iC1<=iClusts;iC1++) {	
      for(iC2=1;iC2<=iClusts;iC2++) {	
	sprintf(strTmp,"%.6f\t",vcInfInter[iC1][iC2]);
	strTab += strTmp;
      }
      strTab += "\n";			
    }
    Write2Log(strTab.c_str());
  } catch(...) {
    Write2Log("Exception in InterClustKLD!!! iC1=%d iC2=%d iClusts=%d",iC1,iC2,iClusts);
    return false;
  }
  return true;
}

prob_t KLDiv(KDTreeHist& p,KDTreeHist& q,bool bAllowZeroQDist)
{	
	int i=0;
#ifdef DO_TIMING
	char sMsg[1024];
	sprintf(sMsg,"KLDiv szp=%d, szq=%d, dims=%d",p.NumElems(),q.NumElems(),p.NumDims());
	ScopedTimer S(sMsg);
#endif
	
	try	//just in case something goes wrong from invalid type of data passed in
	{
		//make sure we can do calculation with div by zero and need same # of dimensions 
		//in each distribution
		if(p.NumElems() < 2 || q.NumElems()<1 || p.NumDims()!=q.NumDims()) 
			return 0.0;

		int iLessP = 0, iLessQ = 0;

		prob_t dpq = 0.0;

		const prob_t eps = 0.0;

		int isz = p.NumElems();

		i = 0;
		for(i=0;i<isz;i++)
		{
			prob_t distp = p.GetNearestRadiusSQ(i,true);
			
			if(distp<=eps)
				continue;

			prob_t distq = q.GetNearestRadiusSQ(p[i],bAllowZeroQDist,true);
			
			if(distq<=eps)
				continue;

			dpq += log2(distq / distp) / 2.0;
		}

		dpq *= ((prob_t)p.NumDims()/p.NumElems());

		dpq += log2( (prob_t)q.NumElems() / (p.NumElems()-1.0 ) );

	#ifdef DO_TIMING
		sprintf(sMsg,"iLessP=%d iLessQ=%d",iLessP,iLessQ);
		MessageBox(0,sMsg,"WClust",MB_ICONINFORMATION);
	#endif
		
		return dpq;
	}
	catch(...)
	{
		char sMsg[1024];
		sprintf(sMsg,"KLDiv caught exception!! i=%d",i);
		fprintf(stderr,"%s\n",sMsg);
		Write2Log(sMsg);
		return -666.0;
	}
}

//resistor avg.
prob_t KLDivSym(KDTreeHist& p,KDTreeHist& q)
{
	prob_t dpq = KLDiv(p,q), dqp = KLDiv(q,p);
	//Write2Log("dpg=%g dqp=%g",dpq,dqp);
	if(!dpq && !dqp) return 0.0;
	return dpq*dqp / (dpq+dqp);
}

prob_t ResistorAvg(prob_t& p,prob_t& q)
{	if(!p && !q) return 0.0f;
	return p*q/(p+q);
}

void CalcUDDist(vector< KDTreeHist >& vDistribs,vector< KDTreeHist >& vCompDistribs,int iClusts,vector<prob_t>& vcInf,vector< vector<prob_t> >& vcInfInter,vector<int>& vCounts)
{
	if(vDistribs.size() != iClusts + 2 || vDistribs.size()!=vCompDistribs.size()) return;

	vcInf = vector<prob_t>(iClusts+1);

	vcInfInter = vector< vector<prob_t> >(iClusts+1);
	int iC=1;
	for(iC=1;iC<=iClusts;iC++) vcInfInter[iC] = vector<prob_t>(iClusts+1);

	//uniqueness from full distribution for each cluster
	int iC1=1,iC2=1;
	for(iC1=1;iC1<=iClusts;iC1++)
		vcInf[iC1] = KLDiv(vDistribs[iC1],vCompDistribs[iC1]);//KLDivSym(vDistribs[iC1],vCompDistribs[iC1]);

	//inter-cluster distinction measures, KL divergence between
	//a cluster and cluster+other_cluster
	for(iC1=1;iC1<=iClusts;iC1++)
	{
		//for(iC2=iC1+1;iC2<=iClusts;iC2++)
		for(iC1=1;iC2<=iClusts;iC2++)
		{
			if(iC1==iC2)
				vcInfInter[iC1][iC2]=0.0;
			else
				//vcInfInter[iC1][iC2] =  vcInfInter[iC2][iC1] = KLDivSym(vDistribs[iC1],vDistribs[iC2]);
				vcInfInter[iC1][iC2] =  KLDiv(vDistribs[iC1],vDistribs[iC2]);
		}
	}
	//add smallest inter-cluster KL-div
	for(iC1=1;iC1<=iClusts;iC1++)
	{
		prob_t dMinInter = 9e10;
		bool bFound = false;
		if(vCounts[iC1])
		for(iC2=1;iC2<=iClusts;iC2++)
		{
			if(iC1 != iC2 && vCounts[iC2] && vcInfInter[iC1][iC2] < dMinInter)
			{
				dMinInter = vcInfInter[iC1][iC2];
				bFound = true;
			}
		}
		if(bFound)
			vcInf[iC1] += dMinInter;
	}
}

double DistSQBG(KDTreeHist& oTree,Neighbor& oN,int iClustID,vector<int>& vClustIDs,int iNodeID)
{
	if(vClustIDs[oN.m_id]!=iClustID)
		return oN.m_dist;
	int iNN = 5;
	double fRad = 5. , fFctr = 5.;
	vector<Neighbor> vnn;
	while(true)
	{
		vnn.resize(iNN);
		oTree.GetKNN(oTree[iNodeID],vnn,iNN,fRad,fFctr);
		int i;
		for(i=0;i<iNN;i++)
			if(vClustIDs[vnn[i].m_id]!=iClustID)
				return vnn[i].m_dist;
		//fRad *= 2.0;
		iNN*=2;
		if(iNN >= oTree.NumElems())
			return 1e5;
	}
}

double DistSQSelf(KDTreeHist& oTree,Neighbor& oN,int iClustID,vector<int>& vClustIDs,int iNodeID)
{
	if(vClustIDs[oN.m_id]==iClustID)
		return oN.m_dist;
	int iNN = 5;
	double fRad = 5. , fFctr = 5.;
	vector<Neighbor> vnn;
	while(true)
	{
		vnn.resize(iNN);
		oTree.GetKNN(oTree[iNodeID],vnn,iNN,fRad,fFctr);
		int i;
		for(i=0;i<iNN;i++)
			if(vClustIDs[vnn[i].m_id]==iClustID)
				return vnn[i].m_dist;
		//fRad *= 2.0;
		iNN*=2;
		if(iNN >= oTree.NumElems())
			return 1e5;
	}
}

double DistD(KDTreeHist& oTree,Neighbor& oN,int iClustID,vector<int>& vClustIDs,int iNodeID)
{
	double d1 = DistSQBG(oTree,oN,iClustID,vClustIDs,iNodeID) ,
		   d2 = DistSQSelf(oTree,oN,iClustID,vClustIDs,iNodeID);
	if(d2) return d1 / d2;
}

prob_t Entropy(KDTreeHist& oTree,vector<Neighbor>& vNeighbors,int iClustID,vector<int>& vClustIDs)
{
	//ScopedTimer S("Entropy...");

	int iFast = 0,iSlow = 0;

	prob_t dEntrop = 0.0;
	int isz = vClustIDs.size() , iV = 0, iT = 0;
	
	double tsz = oTree.NumElems();

	for(iV=0;iV<isz;iV++)
	{
		if(vClustIDs[iV]!=iClustID)
			continue;

		prob_t dDist = 0.0;
		
		Neighbor& oN = vNeighbors[iV];
		if(vClustIDs[oN.m_id]==iClustID && oN.m_dist>0.0)
		{
			dDist = oN.m_dist;
			iFast++;
		}

		if(dDist <= 0.0)
		{
			dDist = oTree.GetNearestRadiusSQ(iT,true);
			iSlow++;
		}

		if(dDist >= 0.0)
		{
			//prob_t dProb = oTree.RProb(MySqrt(dDist));

			//if(dProb >= 0.0)			
			{
				//dEntrop += dProb * log2(dProb);
				dEntrop += log2(tsz*MySqrt(dDist));
			}
		}
		
		iT++;
	}
	//dEntrop *= (oTree.NumDims() / oTree.NumElems());
	dEntrop /= (double) isz;
	//dEntrop += 
	//double r = oTree.NumDims();
	//double N = oTree.NumDims();
	//double PI = 3.14159265;
	//double Sr = r * pow(PI,r/2.0) / Gamma(N/2.0+1.0);
	//double Sr = r * pow(PI,r/2.0) / Gamma(r/2.0+1.0);
	//dEntrop += log2( Sr * (N-1.0) / r);
	return dEntrop;
}

//which 2D slices do we skip?
//slices that are X vs X or
//T1-V(peak) and T1-Peak are highly correlated so skip them too
// when bWC==false, only make sure dimensions aren't identical
inline bool SkipPair(int idx1,int idx2,bool bWC)
{
	if(idx1==idx2) // same dimensions or...
		return true;

	if(!bWC) return false; // no further checks?

	// idx 0 thru 3 == T1,2,3,4-Peak
	// idx 4 thru 7 == T1,2,3,4-V(peak)
	// idx 8 thru 11 == T1,2,3,4-Valley
	// idx 12 thru 15 == T1,2,3,4-V-(valley)

	return (idx1==0 && idx2==4)||   // 0,4 4,0 T1-Peak and T1-V(peak)
		   (idx1==4 && idx2==0)||
		   (idx1==1 && idx2==5)||   // 1,5 5,1 T2-Peak and T2-V(peak)
		   (idx1==5 && idx2==1)||
		   (idx1==2 && idx2==6)||   // 2,6 6,2 T3-Peak and T3-V(peak)
		   (idx1==6 && idx2==2)||		   
		   (idx1==3 && idx2==7)||   // 3,7 7,3 T4-Peak and T4-V(peak)
		   (idx1==7 && idx2==3)||
		   (idx1==8 && idx2==12)||  // 8,12 12,8 T1-Valley and T1-V(valley)
		   (idx1==12 && idx2==8)||
		   (idx1==9 && idx2==13)||  // 9,13 13,9 T2-Valley and T2-V(valley)
		   (idx1==13 && idx2==9)||
		   (idx1==10 && idx2==14)|| // 10,14 14,10 T3-Valley and T3-V(valley)
		   (idx1==14 && idx2==10)||
		   (idx1==11 && idx2==15)|| // 11,15 15,11 T4-Valley and T4-V(valley)
		   (idx1==15 && idx2==11);	
}
// HasSkipPair - check if ints in s and i form a "SkipPair"
bool HasSkipPair(set<int>& s,int i,bool bWC)
{
  set<int>::iterator IT = s.begin();
  for(;IT!=s.end();IT++)
    if(SkipPair(*IT,i,bWC))
      return true;
  return false;
}

float KLDivSym1D(vector<float>& vFloat,vector<int>& vClustIDs,vector<int>& vCounts,int iCols,int iRows,int iClust,int iDim)
{
	vector<float> v1DFloatClust(vCounts[iClust]),v1DFloatComp(iRows-vCounts[iClust]);
	int idxClust = 0, idxComp = 0;
	KDTreeHist o1DTClust,o1DTComp;
	int iV = 0;
	for(iV=0;iV<vClustIDs.size();iV++)
	{	if(vClustIDs[iV]==iClust)
			v1DFloatClust[idxClust++]=vFloat[iV*iCols+iDim];
		else
			v1DFloatComp[idxComp++]=vFloat[iV*iCols+iDim];
	}
	o1DTClust.SetData(1,&v1DFloatClust[0],vCounts[iClust]);
	o1DTComp.SetData(1,&v1DFloatComp[0],iRows-vCounts[iClust]);
	return KLDivSym(o1DTClust,o1DTComp);
}

bool FindBest2DDims(vector<float>& vFloat,vector<float>& vRange,int iClusts,int iCols,int iBestDims,vector<int>& vCounts,vector<int>& vClustIDs,A2D<int>& vBestDims,A2D<KLD2D>& vKLDivs,bool bWC,bool bTime,vector<ClusterInfo>& vInfo)
{	vBestDims.Init(iClusts+1,iBestDims);//each cluster will get iBestDims to perform multidimensional kldiv on later
	bool bInit = false;
	//vKLDivs.Init(iClusts+1,iBestDims);
	int iC = 1 , iRows = vClustIDs.size();
	double dJnk = 0.0;
	const float fMinRange = 0.009; //min range for a dimension to be usable
	const float fMaxRange = 1e7; //max range for a dimension to be usable
	char msg[2048];
	CStopWatch oTimer;
	for(iC=1;iC<=iClusts;iC++)
        {	if(bTime) oTimer.startTimer(); // timer for finding best 2D slices
		vector<KLD2D> vKLDivTmp(IntegerSum(iCols-1));
		int iD1,iD2, iK = 0;
		for(iD1=0;iD1<iCols;iD1++)
		{	
			for(iD2=iD1+1;iD2<iCols;iD2++,dJnk++)
			{
			        if(SkipPair(iD1,iD2,bWC)) continue;

				//exclude 2D slice consisting of empty signal -- occurs when a wire is grounded
				//also exclude dimensions where the range is so huge that its likely to be noise, this rarely
				//happens but when it does it can produce spurious results by forcing together very tightly points
				//that shouldn't be so close
				if(vRange[iD1]<fMinRange || vRange[iD2]<fMinRange || vRange[iD1]>fMaxRange || vRange[iD2]>fMaxRange)
				{	Write2Log("Skipping slice %d %d with ranges %.12f %.12f",iD1,iD2,vRange[iD1],vRange[iD2]);
					continue;
				}				
				vector<float> v2DFloatClust(2*vCounts[iC]),v2DFloatComp(2*(iRows-vCounts[iC]));
				int idxClust = 0, idxComp = 0;
				KDTreeHist o2DTClust,o2DTComp;
				int iV = 0;
				for(iV=0;iV<vClustIDs.size();iV++)
				{	//initialize the trees
					if(vClustIDs[iV]==iC)
					{	v2DFloatClust[idxClust++]=vFloat[iV*iCols+iD1];
						v2DFloatClust[idxClust++]=vFloat[iV*iCols+iD2];
					}
					else
					{	v2DFloatComp[idxComp++]=vFloat[iV*iCols+iD1];
						v2DFloatComp[idxComp++]=vFloat[iV*iCols+iD2];
					}
				}
				o2DTClust.SetData(2,&v2DFloatClust[0],vCounts[iC]);
				o2DTComp.SetData(2,&v2DFloatComp[0],iRows-vCounts[iC]);
				prob_t kld = KLDivSym(o2DTClust,o2DTComp); //compute kld
				//Write2Log("2D kld C%d dims:%s %s = %.4f",iC,*vAxes[iD1],*vAxes[iD2],kld);
				vKLDivTmp[iK++].Init(iD1,iD2,kld);//store kldiv
			}
		}
		vKLDivTmp.resize(iK);
		//sort results by kldiv values
		sort(vKLDivTmp.begin(),vKLDivTmp.end());

		if(!bInit)
		{
			vKLDivs.Init(iClusts+1,iK); // init here
			bInit=true;
		}

		copy(vKLDivTmp.begin(),vKLDivTmp.end(),&vKLDivs[iC][0]); // copy for later use

		//go through results, picking out top 8 dimensions, make sure not to pick
		//the same dimension twice
		int iFound = 0;
		set<int> sDims; // stores dimensions already picked
		int idx = iK-1; // start at best kldiv pair
		set< pair<float,int> > vKD1D;
		for(;idx>=0 && iFound<iBestDims;idx--)
		{	
   		           bool bHas1=sDims.find(vKLDivTmp[idx].m_iD1)!=sDims.end() || HasSkipPair(sDims,vKLDivTmp[idx].m_iD1,bWC),
  			        bHas2=sDims.find(vKLDivTmp[idx].m_iD2)!=sDims.end() || HasSkipPair(sDims,vKLDivTmp[idx].m_iD2,bWC);
			if(iFound+2<=iBestDims)
			{	//can add both
				if(!bHas1)
				{	sDims.insert(vKLDivTmp[idx].m_iD1);	//store which dims we already have
					vBestDims[iC][iBestDims-iFound-1]=vKLDivTmp[idx].m_iD1;
					iFound++;
				}
				if(!bHas2)
				{	sDims.insert(vKLDivTmp[idx].m_iD2);
					vBestDims[iC][iBestDims-iFound-1]=vKLDivTmp[idx].m_iD2;
					iFound++;
				}
			}
			else if(iFound+1<=iBestDims)
			{	//can add only 1
				if(bHas1 && !bHas2)
				{	sDims.insert(vKLDivTmp[idx].m_iD2);
					vBestDims[iC][iBestDims-iFound-1]=vKLDivTmp[idx].m_iD2;
					iFound++;
				}
				else if(bHas2 && !bHas1)
				{	sDims.insert(vKLDivTmp[idx].m_iD1);	//store which dims we already have
					vBestDims[iC][iBestDims-iFound-1]=vKLDivTmp[idx].m_iD1;
					iFound++;
				}
				else if(!bHas1)
				{	sDims.insert(vKLDivTmp[idx].m_iD1);	//store which dims we already have
					vBestDims[iC][iBestDims-iFound-1]=vKLDivTmp[idx].m_iD1;
					iFound++;
				}
				else if(!bHas2) // can't ever go here
				{	sDims.insert(vKLDivTmp[idx].m_iD2);
					vBestDims[iC][iBestDims-iFound-1]=vKLDivTmp[idx].m_iD2;
					iFound++;
				}
				/*               fun with truth tables
						has1   has2  !has1  !has2   has1 && !has2   has2 && !has1
						true   true  false  false       false          false
						true  false  false  true        true           false
						false  true  true   false       false          true
						false false  true   true        false          false
				*/
			}
		}
#if 0
		Write2Log("\nClust%d 2D kldiv pairs(best 16) info follows:\n",iC);
		LogF F;
		FILE* fp = F.Open();
		int y=iK-16>=0?iK-16:0;
		for(;y<iK;y++)
		{	fprintf(fp,"pair%d D1=%d D2=%d kld=%.4f\n",
				y,vKLDivTmp[y].m_iD1,vKLDivTmp[y].m_iD2,vKLDivTmp[y].m_kld);
		} fprintf(fp,"\n\n");
#endif
		if(bTime) { // save time to find best 2D slices for this cluster
		  oTimer.stopTimer();
		  vInfo[iC].m_dDimTime = oTimer.getElapsedTime();
		}
	}
	return true;
}

prob_t MaxAbsCorr(vector<vector<prob_t> >& vCorrelMat,int* vBestDims,int iDim,int iBestDims)
{
	int i;
	prob_t bestCor = 0.0, tmpCor = 0.0;
	for(i=iBestDims-1;i>=0;i--)
		if(vBestDims[i]!=iDim && (tmpCor=fabs(vCorrelMat[vBestDims[i]][iDim]))>bestCor)
			bestCor=tmpCor;
	return bestCor;
}

prob_t MaxAbsCorr(vector<vector<prob_t> >& vCorrelMat,vector< pair<float,int> >& vKLDDim,int iDim)
{
	int iSz = vKLDDim.size(), i = 0;
	prob_t maxScore = 0.0 , tmpScore = 0.0;
	for(i=0;i<iSz;i++)
	{	if( vKLDDim[i].second == iDim ) continue;
		if( (tmpScore = fabs(vCorrelMat[vKLDDim[i].second][iDim]))>maxScore)
			maxScore = tmpScore;
	}
	return maxScore;
}

prob_t KLDCorVal(vector<vector<prob_t> >& vCorrelMat,vector< pair<float,int> >& vKLDDim)
{
	int iSz = vKLDDim.size() , i = 0;
	prob_t score = 0.0;
	for(i=0;i<iSz;i++)
		score += vKLDDim[i].first * (1.0 - MaxAbsCorr(vCorrelMat,vKLDDim,vKLDDim[i].second));
	return score;
}

int GetReplaceIndex(vector<vector<prob_t> >& vCorrelMat,vector< pair<float,int> >& vOrig,pair<float,int>& oP)
{
	int iSz = vOrig.size() , iBestIDX = -1 , i = 0;
	prob_t maxVal = KLDCorVal(vCorrelMat,vOrig);
	prob_t curVal = maxVal;
	for(i=0;i<iSz;i++)
	{
		vector< pair<float,int> > vTmp(vOrig);
		vTmp[i]=oP;
		if( (curVal=KLDCorVal(vCorrelMat,vTmp)) > maxVal )
		{
			maxVal = curVal;
			iBestIDX = i;
		}
	}
	return iBestIDX;
}

bool FindBest1DDims(vector<float>& vFloat,int iClusts,int iCols,int iBestDims,vector<int>& vCounts,vector<int>& vClustIDs,A2D<int>& vBestDims,A2D<prob_t>& vKLDivs)
{	vBestDims.Init(iClusts+1,iBestDims);
	vKLDivs.Init(iClusts+1,iBestDims);
	vKLDivs.Fill(-99999.9f);
	int iC = 1 , iRows = vClustIDs.size();
	double dJnk = 0.0;
	char msg[2048];
	vector< vector<float> > vCorrelMat;
	vector<float> vMean;
	sprintf(msg,"Computing %d X %d correlation matrix",iCols,iCols);
	fprintf(stderr,"%s\n",msg);
	CovarMat(vFloat,vClustIDs.size(),iCols,vCorrelMat,vMean,true);
	for(iC=1;iC<=iClusts;iC++)
	{	int iD;
		vector< pair<float, int> > vKLDDims(iCols);
		for(iD=0;iD<iCols;iD++,dJnk++)
		  {	sprintf(msg,"Finding best %d dimensions for cluster %d of %d : Dim=%d",iBestDims,iC,iClusts,iD);
		        fprintf(stderr,"%s\n",msg);
			vector<float> v1DFloatClust(vCounts[iC]),v1DFloatComp(iRows-vCounts[iC]);
			int idxClust = 0, idxComp = 0;
			KDTreeHist o1DTClust,o1DTComp;
			int iV = 0;
			for(iV=0;iV<vClustIDs.size();iV++)
			{	if(vClustIDs[iV]==iC)
					v1DFloatClust[idxClust++]=vFloat[iV*iCols+iD];
				else
					v1DFloatComp[idxComp++]=vFloat[iV*iCols+iD];
			}
			o1DTClust.SetData(1,&v1DFloatClust[0],vCounts[iC]);
			o1DTComp.SetData(1,&v1DFloatComp[0],iRows-vCounts[iC]);
			prob_t kld = KLDivSym(o1DTClust,o1DTComp);
			
			vKLDDims[iD]=pair<float,int>(kld,iD);
			/*pair<float,int> oP(kld,iD);
			if(vKLDDims.size()<iBestDims)
				vKLDDims.push_back(oP);
			else
			{	int idx = GetReplaceIndex(vCorrelMat,vKLDDims,oP);
				if(idx != -1)
					vKLDDims[idx]=oP;
			}*/
		}
		/*sort(vKLDTmp.begin(),vKLDTmp.end());
		int iJ , idx = iBestDims-1, iFound=0;
		for(iJ=iCols-1;iJ>=0;iJ--)
		{
			int iK;
			bool bRedund = false;
			for(iK=iBestDims-1;iK>=idx;iK--)
			{
				if(SkipPair(vBestDims[iC][iK],vKLDTmp[iJ].second))
				{
					bRedund = true;
					break;
				}
			}
			if(!bRedund)
			{
				vBestDims[iC][idx]=vKLDTmp[iJ].second;
				vKLDivs[iC][idx]=vKLDTmp[iJ].first;
				idx--;
				iFound++;
			}
			if(iFound==iBestDims)
				break;
		}*/
		vector<int> dimIDs(iCols) , bestdimIDs(iCols);
		int iJ = 0;
		for(iJ=0;iJ<iCols;iJ++) dimIDs[iJ]=iJ;
		//LogF F; FILE* fp=F.Open();
		float maxScore = -10000.0 , tmpScore = 0.0;
		btb::combination_init(&dimIDs[0],&dimIDs[iBestDims],&dimIDs[dimIDs.size()]);
		do
		{	vector< pair<float,int> > vKLDTmp(iBestDims);
			int iK = 0;
			for(;iK<iBestDims;iK++)
			{
				vKLDTmp[iK]=vKLDDims[dimIDs[iK]];
				//fprintf(fp,"%d\t",dimIDs[iK]);
			} //fprintf(fp,"\n");
			if((tmpScore=KLDCorVal(vCorrelMat,vKLDTmp))>maxScore)
			{   maxScore = tmpScore;
				bestdimIDs=dimIDs;
				Write2Log("maxScore=%.2f",maxScore);
				// WriteVec2Log(bestdimIDs);
			}

		}while(btb::next_combination(&dimIDs[0], &dimIDs[iBestDims], &dimIDs[dimIDs.size()]));
		//F.Close();
		vector< pair<float,int> > vKLDTmp(iBestDims);
		for(iJ=0;iJ<iBestDims;iJ++) vKLDTmp[iJ] = vKLDDims[bestdimIDs[iJ]];
		vKLDDims=vKLDTmp;
		sort(vKLDDims.begin(),vKLDDims.end());
		for(iJ=0;iJ<vKLDDims.size();iJ++)
		{	vBestDims[iC][iJ]=vKLDDims[iJ].second;
			vKLDivs[iC][iJ]=vKLDDims[iJ].first;
		}
	}
	return true;
}

#ifdef _OLDER_V
bool FindBest1DDims(vector<float>& vFloat,int iClusts,int iCols,int iBestDims,vector<int>& vCounts,vector<int>& vClustIDs,A2D<int>& vBestDims,A2D<prob_t>& vKLDivs)
{	vBestDims.Init(iClusts+1,iBestDims);
	vKLDivs.Init(iClusts+1,iBestDims);
	vKLDivs.Fill(-99999.9f);
	int iC = 1 , iRows = vClustIDs.size();
	double dJnk = 0.0;
	char msg[2048];
	for(iC=1;iC<=iClusts;iC++)
	{	int iD;
		for(iD=0;iD<iCols;iD++,dJnk++)
		  {	sprintf(msg,"Finding best %d dimensions for cluster %d of %d : Dim=%d",iBestDims,iC,iClusts,iD);
		        fprintf(stderr,"%s\n",msg);
			vector<float> v1DFloatClust(vCounts[iC]),v1DFloatComp(iRows-vCounts[iC]);
			int idxClust = 0, idxComp = 0;
			KDTreeHist o1DTClust,o1DTComp;
			int iV = 0;
			for(iV=0;iV<vClustIDs.size();iV++)
			{	if(vClustIDs[iV]==iC)
					v1DFloatClust[idxClust++]=vFloat[iV*iCols+iD];
				else
					v1DFloatComp[idxComp++]=vFloat[iV*iCols+iD];
			}
			o1DTClust.SetData(1,&v1DFloatClust[0],vCounts[iC]);
			o1DTComp.SetData(1,&v1DFloatComp[0],iRows-vCounts[iC]);
			prob_t kld = KLDivSym(o1DTClust,o1DTComp);
			int iJ;
			for(iJ=iBestDims-1;iJ>=0;iJ--)
			{	if(kld>vKLDivs[iC][iJ])
				{	int iJ2;
					for(iJ2=0;iJ2<iJ;iJ2++)
					{	vKLDivs[iC][iJ2]=vKLDivs[iC][iJ2+1];
						vBestDims[iC][iJ2]=vBestDims[iC][iJ2+1];
					}
					vKLDivs[iC][iJ]=kld;
					vBestDims[iC][iJ]=iD;
					break;
				}
			}
		}
	}
	return true;
}
#endif
