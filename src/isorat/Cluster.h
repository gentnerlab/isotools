//////////////////////////////////////////////////////////////////////
// Cluster.h: interface for the CCluster class.
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

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include "MyObj.h"
#include "A2D.h"

struct ClusterInfo {
  int m_iMyID;
  int m_iSz;	//# of spikes loaded in
  prob_t m_dIsolationDist;
  prob_t m_dLRatio;
  double m_dIsoDTime;
  double m_dLRatTime;
  
  ClusterInfo()
  :m_iMyID(0),
   m_iSz(0),
   m_dIsolationDist(-1.0), // starts invalid
   m_dLRatio(-1.0), // starts invalid
   m_dIsoDTime(0.0),
   m_dLRatTime(0.0)
  {}
	
};

////////////////////////////////////////////////////////////////////////
// CCluster
class CCluster : public CMyObject {
public:	
  vector< ClusterInfo > m_vInfo;   //cluster info. 
  int m_iNumClusts;  
  void GetClusterInfo();
public:
  CCluster();
  virtual ~CCluster() {}
  void	Clear();
  int		GetNumClusts();
};

#endif 
