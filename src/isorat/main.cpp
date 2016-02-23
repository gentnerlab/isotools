// $Id: main.cpp,v 1.1 2011/03/07 16:36:01 samn Exp $ 

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
//////////////////////////////////////////////////////////////////////

#include "Vertex.h"

#include <stdio.h>
#include <string.h>

bool ClustIsolationDist(CCluster& Clusters,CVerxStack& oS,vector<double>& vFloat,vector<int>& vClustIDs,vector<int>& vClustCounts,int iClusts,int iRows,int iCols,bool bTime);
bool ClustLRatio(CCluster& Clusters,CVerxStack& oS,vector<double>& vFloat,vector<int>& vClustIDs,vector<int>& vClustCounts,int iClusts,int iRows,int iCols,bool bTime);


void prhelp() {

  fprintf(stderr,"\nisorat: copyright (C) 2003-2011 Sam Neymotin & BioSignal Group.\n\n");

  fprintf(stderr,"This program runs the IsolationDistance and LRatio measures on cluster feature files.\n\n");

  fprintf(stderr,"Command line usage: isorat input output [-time]\n\n");

  fprintf(stderr," The input file is in ASCII/text format.\n");
  fprintf(stderr,"  The first line of input is a header for the names of the dimensions on the following lines.\n");
  fprintf(stderr,"  The dimension names are separated by white-space.\n");
  fprintf(stderr,"  Note that the first dimension must be the cluster identifier.\n");
  fprintf(stderr,"  After the header, the remaining lines of the input file are feature vector records.\n");
  fprintf(stderr,"  Each feature vector record consists of: clusterid as an integer followed by decimal values\n");
  fprintf(stderr,"  of the feature vectors, separated by spaces. Each record is on a separate line.\n");
  fprintf(stderr,"  NOTE that a clusterid of 0 indicates the noise distribution. Clusters that correspond\n");
  fprintf(stderr,"  to isolated single units (non-noise) should have a clusterid greater or equal to 1.\n\n");

  fprintf(stderr," The output is saved as a text file with one line for each input cluster.\n");
  fprintf(stderr,"  Each line consists of: IsoD (cluster separation from background)\n");
  fprintf(stderr,"  and LRatio (cluster separation from background).\n\n");

  fprintf(stderr," The -time flag specifies to perform timing test, if set, will add time column for each cluster.\n\n");
}

int main (int argc, char** argv) {
  if(argc < 3) {prhelp();  return 0;}
  CVerxStack cs;  
  if(cs.ReadAsciiData(argv[1])) 
  	printf("read %d %d-dim vectors from %s\n",cs.m_NumVerx,cs.m_iNumDim,argv[1]);
  else {
  	fprintf(stderr,"couldn't read vectors from %s!\n",argv[1]);
  	return 0;
  }  
  bool bTime = argc >= 4 && !strcmp(argv[3],"-time") ? true : false; // whether to do timing test
  CCluster& oC = cs.m_oClusters;
  vector<double> vDat; vector<float> vRange;
  int iClusts=oC.GetNumClusts();
  vector<int> vCounts(iClusts+1);
  cs.GetCounts(vCounts,iClusts);	

  //multidimensional distributions - continuous
  int iRows=0,iCols=0;
  cs.GetDataV(iRows,iCols,vDat,vRange);

  vector<int> vClustIDs;
  cs.GetClustIDs(vClustIDs);	

  if(!ClustIsolationDist(oC,cs,vDat,vClustIDs,vCounts,iClusts,iRows,iCols,bTime) ||
     !ClustLRatio(oC,cs,vDat,vClustIDs,vCounts,iClusts,iRows,iCols,bTime)) {
    fprintf(stderr,"isorat ERR: couldn't calculate IsoD/LRatio!\n");
    return 0;
  }
  int iRet = 1;
  FILE* fp = fopen(argv[2],"w");
  if(!fp) {
  	fprintf(stderr,"couldn't open output file %s!\n",argv[2]);
  	iRet = 0;
  }
  printf("Cluster info:\n");
  int i;
  if(bTime) {
    for(i=1;i<=iClusts;i++) {
      ClusterInfo& oCI = oC.m_vInfo[i];
      printf("\tC%d: NumSpikes=%d, IsoD=%g, LRatio=%g, IsoDTime=%g, LRatTime=%g\n",
	     i,vCounts[i],oCI.m_dIsolationDist,oCI.m_dLRatio,oCI.m_dIsoDTime,oCI.m_dLRatTime);
      if(fp) fprintf(fp,"%g %g %g %g\n",oCI.m_dIsolationDist,oCI.m_dLRatio,oCI.m_dIsoDTime,oCI.m_dLRatTime);
    }
  } else for(i=1;i<=iClusts;i++) {
  	ClusterInfo& oCI = oC.m_vInfo[i];
  	printf("\tC%d: NumSpikes=%d, IsoD=%g, LRatio=%g\n",
	       i,vCounts[i],oCI.m_dIsolationDist,oCI.m_dLRatio);
  	if(fp) fprintf(fp,"%g %g\n",oCI.m_dIsolationDist,oCI.m_dLRatio);
  }
  if(fp) fclose(fp);

  return iRet;
}
