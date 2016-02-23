// $Id: main.cpp,v 1.7 2011/08/01 13:56:01 samn Exp $ 
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

#include "KDTree.h"
#include "Vertex.h"
#include "InfoT.h"
#include "hr_time.h"

#include <stdio.h>

bool CalcClustInfo(CVerxStack& DataStack,bool bSymmetric,A2D<Neighbor>& vnn,
		   vector<float>& vFloat,KDTreeHist& oTree,vector<float>& vRange,bool bWC,bool bTime);

void prhelp() {

  fprintf(stderr,"\nisoi: copyright (C) 2003 - 2011, Sam Neymotin and BioSignalGroup.\n\n");
  fprintf(stderr,"This program runs the Isolation Information measure on cluster feature files.\n\n");

  fprintf(stderr,"Command-line usage: isoi input output [-wc][-time]\n\n");

  fprintf(stderr,"The input file is in ASCII/text format.\n");
  fprintf(stderr,"The first line of input is a header for the names of the dimensions on the following lines.\n");
  fprintf(stderr,"The dimension names are separated by white-space.\n");
  fprintf(stderr,"Note that the first dimension must be the cluster identifier.\n");
  fprintf(stderr,"After the header, the remaining lines of the input file are feature vector records.\n");
  fprintf(stderr,"Each feature vector record consists of: clusterid as an integer followed by decimal values\n");
  fprintf(stderr,"of the feature vectors, separated by spaces. Each record is on a separate line.\n");
  fprintf(stderr,"NOTE that a clusterid of 0 indicates the noise distribution. Clusters that correspond\n");
  fprintf(stderr,"to isolated single units (non-noise) should have a clusterid greater or equal to 1.\n\n");

  fprintf(stderr,"The output is saved as a text file with one line for each input cluster.\n");
  fprintf(stderr,"Each line consists of: IsoI_BG (cluster separation from background),IsoI_NN (cluster separation\n");
  fprintf(stderr,"from nearest neighboring cluster), ClosestCluster (nearest neighboring cluster identifier).\n\n");

  fprintf(stderr,"The -wc flag specifies that input features are the same as those used in the WClust spike-sorting\n");
  fprintf(stderr,"program. When this flag is provided, 2D slices consisting of highly correlated dimensions, such as\n");
  fprintf(stderr,"T1-Peak and T1-V(peak) are discarded. If you are not using the -wc flag, it's up to you to ensure\n");
  fprintf(stderr,"the dimensions provided to isoi are not overly correlated.\n\n");

  fprintf(stderr,"The -time flag specifies to perform timing test, if set, will add time column for each cluster.\n\n");
}

int main (int argc, char** argv) {
  if(argc < 3) {prhelp();  return 0;}
  bool bWC = false, bTime = false;
  if(argc >= 4 && !strcmp(argv[3],"-wc")) bWC=true;
  if(argc >= 5 && !strcmp(argv[4],"-time")) bTime=true;
  CVerxStack cs;  
  if(cs.ReadAsciiData(argv[1])) 
    printf("read %d %d-dim vectors from %s\n",cs.m_NumVerx,cs.m_iNumDim,argv[1]);
  else {
    fprintf(stderr,"couldn't read vectors from %s!\n",argv[1]);
    return 0;
  }  
  CCluster& oC = cs.m_oClusters;
  A2D<Neighbor> vnn;
  vector<float> vFloat,vRange;
  KDTreeHist oTree;
  if(!CalcClustInfo(cs,true,vnn,vFloat,oTree,vRange,bWC,bTime)) {
    fprintf(stderr,"couldn't calculate cluster info!\n");
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
    for(i=1;i<=oC.GetNumClusts();i++) {
      ClusterInfo& oCI = oC.m_vInfo[i];
      printf("\tC%d: NumSpikes=%d, IsoI_BG=%g, IsoI_NN=%g, ClosestClust=%d, DimTime=%g, BGTime=%g, NNTime=%g\n",
	     i,oCI.m_iSz,oCI.m_fBGInfoGain,oCI.m_fInterClustGain,oCI.m_iClosestID,oCI.m_dDimTime,oCI.m_dBGTime,oCI.m_dNNTime);
      if(fp) fprintf(fp,"%g %g %d %g %g %g\n",oCI.m_fBGInfoGain,oCI.m_fInterClustGain,oCI.m_iClosestID,
		     oCI.m_dDimTime,oCI.m_dBGTime,oCI.m_dNNTime);
    }
  } else for(i=1;i<=oC.GetNumClusts();i++) {
  	ClusterInfo& oCI = oC.m_vInfo[i];
  	printf("\tC%d: NumSpikes=%d, IsoI_BG=%g, IsoI_NN=%g, ClosestClust=%d\n",
	       i,oCI.m_iSz,oCI.m_fBGInfoGain,oCI.m_fInterClustGain,oCI.m_iClosestID);
  	if(fp) fprintf(fp,"%g %g %d\n",oCI.m_fBGInfoGain,oCI.m_fInterClustGain,oCI.m_iClosestID);
  }
  if(fp) fclose(fp);
  return iRet;
}
