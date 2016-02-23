// $Id: Cluster.h,v 1.6 2011/08/01 13:55:47 samn Exp $ 
// Cluster.h: interface for the CCluster class.
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

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include "MyObj.h"
#include "A2D.h"

struct KLD2D {
  int m_iD1;
  int m_iD2;
  prob_t m_kld;
  KLD2D() {
    m_iD1=m_iD2=-1;
    m_kld=-66666.6f;
  }
  void Init(int iD1,int iD2,prob_t kld) {
    m_iD1=iD1;
    m_iD2=iD2;
    m_kld=kld;
  }
};

inline bool operator<(KLD2D& k1,KLD2D& k2) {
  return k1.m_kld<k2.m_kld;
}

inline bool operator<(const KLD2D& k1,const KLD2D& k2) {
  return k1.m_kld<k2.m_kld;
}

struct ClusterInfo
{
  prob_t m_fBGInfoGain;
  prob_t m_fInterClustGain;
  int m_iMyID;
  int m_iClosestID;
  prob_t m_fPrctKNNInClust;
  int m_iSz;	//# of spikes loaded in
  double m_dDimTime; // time to find best dims
  double m_dBGTime; // time to get IsoI_BG
  double m_dNNTime; // time to get IsoI_NN
	
  ClusterInfo()
    :m_fBGInfoGain(0.0),
     m_fInterClustGain(0.0),
     m_iMyID(0),
     m_iClosestID(0),
     m_fPrctKNNInClust(0.0),
     m_iSz(0),
     m_dDimTime(0.0),
     m_dBGTime(0.0),
     m_dNNTime(0.0)
     {}
	
  ClusterInfo(prob_t fBGInfoGain,prob_t fInterClustGain,int iID,int iClosestID)
  :m_fBGInfoGain(fBGInfoGain),
   m_fInterClustGain(fInterClustGain),
   m_iMyID(iID),
   m_iClosestID(iClosestID),
   m_fPrctKNNInClust(0.0),
   m_iSz(0),
   m_dDimTime(0.0),
   m_dBGTime(0.0),
   m_dNNTime(0.0)
   {}
};

//struct storing cluster quality options
struct CQOpts {
  int		m_iNNToFind; //how many nearest neighbors to search for initially (when computing KLDIV)
  bool		m_bFastKLD;  //whether to use fast KLDIV (uses more mem.)
  bool		m_bInterKLD; //whether to compute inter-cluster kldiv
  bool		m_bDoKLDiv;  //whether to compute kldiv at all
  int		m_iBestDims;	//how many dimensions to use for kldiv
  bool		m_bFindBestDims; //whether to find best dimensions before computing KLDIV
  bool          m_b1DBest;		//iff==true use 1D KLDivs, otherwise use 2D KLDivs for finding best

  void Default() {
    m_iNNToFind = 1;
    m_bFastKLD = false;
    m_bInterKLD = true;
    m_bDoKLDiv = true;
    m_iBestDims = 8;
    m_bFindBestDims = true;
    m_b1DBest = false; //use 2D by default
  }
  CQOpts(){ Default(); }
};

////////////////////////////////////////////////////////////////////////
// CCluster
class CCluster : public CMyObject  {
public:	
  CQOpts m_oCQO;
  A2D<int>	m_vBestDims;	//best dimensions for each cluster
  A2D<prob_t>	m_v1DKLDivs;	//1D kldivs (clust vs background) of best dimensions for each cluster
  A2D<KLD2D>  m_v2DKLDivs;    //2D kldivs (clust vs background) of all pairs of dimensions
  vector< ClusterInfo > m_vInfo;   //cluster info. 
  int m_iNumClusts;
  void GetClusterInfo();
public:
  CCluster();
  virtual ~CCluster() {}
  void	Clear();
  int	GetNumClusts();
};

#endif 
