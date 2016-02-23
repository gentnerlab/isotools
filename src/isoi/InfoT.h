// $Id: InfoT.h,v 1.5 2011/08/01 13:56:10 samn Exp $ 
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

#ifndef INFOT_H
#define INFOT_H

#include "Hist.h"

//calculates kldiv from p to q using previously computed nearest neighbors
//this version assumes both p and q have ids of iClustIDP and iClustIDQ
//this means q is not the 'background' cluster but an actual cluster!!
prob_t FastKLDivPQ(KDTreeHist& p,KDTreeHist& q,Neighbor** vNeighbors,int iClustIDP,int iClustIDQ,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount);

//calculates kldiv from p to notp using previously computed nearest neighbors
//this version assumes p is an actual cluster and notp is
//the complement of p == 'background' cluster -- doesn't have it's own ID!!
prob_t FastKLDivPNOTP(KDTreeHist& p,KDTreeHist& notp,Neighbor** vNeighbors,int iClustIDP,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount);

//calculates kldiv from notp to p using previously computed nearest neighbors
//notp is the complement of p -- or the background cluster relative to p -- doesn't have it's own ID!
prob_t FastKLDivNOTPP(KDTreeHist& p,KDTreeHist& notp,Neighbor** vNeighbors,int iClustIDP,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount);

//calculates (symmetric) resistor average of kldiv. (kldiv(p,q)*kldiv(q,p))/(kldiv(p,q)+kldiv(q,p))
//here both p and q are actual clusters with proper IDs.
//neither p nor q is a 'background' cluster relative to the other
prob_t FastKLDivSymPQ(KDTreeHist& p,KDTreeHist& q,Neighbor** vNeighbors,int iClustIDP,int iClustIDQ,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount);

//calculates (symmetric) resistor average of kldiv between p and not p
//here p is an actual cluster with proper ID and not p is the background cluster relative to p
prob_t FastKLDivSymPNOTP(KDTreeHist& p,KDTreeHist& notp,Neighbor** pNeighbors,int iClustIDP,vector<int>& vClustIDs,int iNNToFind,vector<int>& vNCount);

//thread to calculate minimum inter-cluster kldiv for each cluster and add it to Clusters.m_vInfo
bool InterClustKLD(CCluster& Clusters,vector<KDTreeHist>& vDistribs,vector<int>& vClustIDs,vector<int>& vClustCounts,int iClusts,bool bFast,Neighbor** vnn,int WhichDraw,int iNNToFind,vector<int>& vNCount,vector<float>& vFloat,int iRows,int iCols,bool bTime);

//get prct of iNNToFind neighbors that have same cluster for each cluster
bool NeighborClustCount(Neighbor** vNeighbors,vector<int>& vClustIDs,vector<ClusterInfo>& vPrct,int iNNToFind,int iClusts);

//kldiv
//bAllowZeroQDist is used for when you don't want to allow any zero distances for
//probability calculations of elements in distribution q. this is not generally useful
//for distributions containing different elements, but if q is the FULL background distribution
//bAllowZeroQDist should be set to false. default == true and is whats used for standard KLDiv
//calculations of cluster quality
prob_t KLDiv(KDTreeHist& p,KDTreeHist& q,bool bAllowZeroQDist=true);
//symmetric kldiv
prob_t KLDivSym(KDTreeHist& p,KDTreeHist& q);
//resistor average
prob_t ResistorAvg(prob_t& p,prob_t& q);

//entropy using previously computed nearest neighbors
prob_t Entropy(KDTreeHist& oTree,vector< Neighbor >& vNeighbors,int iClustID,vector<int>& vClustIDs);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                FindBestDims                                                          //
//                 find best dimensions for each cluster using KLDIV of 1 or 2 dimensional slices and ranking           //
//                                                                                                                      //
//  vFloat is full data matrix                                                                                          //
//  iClusts == # of clusters                                                                                            //
//  iCols == full # of dimensions                                                                                       //
//  iBestDims == desired # of dimensions after selection                                                                //
//  vCounts == vector of # of vectors in a cluster                                                                      //
//  vClustIDs == cluster IDs for each vector                                                                            //
//  vBestDims == stores best dimensions for all clusters                                                                //
//  vKLDivs == stores kldiv values when using 1D slices                                                                 //
//  b1D == true iff user wants 1D slices to estimate best dimensions                                                    //
//                                                                                                                      //
bool FindBestDims(vector<float>& vFloat,int iClusts,int iCols,int iBestDims,vector<int>& vCounts,                       //
		  vector<int>& vClustIDs,A2D<int>& vBestDims,A2D<prob_t>& vKLDivs,bool b1D);
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FindBest2DDims - bWC argument is true when using tetrode and same dimensions/order of dimensions as WClust
// when bWC is true, checks for redundant dimensions, i.e. T1-Peak and T1-V(peak). If bWC is not true, it's up to
// the user to ensure dimensions are not overly redundant. Can put in threshold for correlation btwn dimensions in
// a 2D slice later...
bool FindBest2DDims(vector<float>& vFloat,vector<float>& vRange,int iClusts,int iCols,int iBestDims,vector<int>& vCounts,
		    vector<int>& vClustIDs,A2D<int>& vBestDims,A2D<KLD2D>& vKLDivs,bool bWC,bool bTime,vector<ClusterInfo>& vInfo);
// FindBest1DDims
bool FindBest1DDims(vector<float>& vFloat,int iClusts,int iCols,int iBestDims,vector<int>& vCounts,
		    vector<int>& vClustIDs,A2D<int>& vBestDims,A2D<prob_t>& vKLDivs);

#endif
