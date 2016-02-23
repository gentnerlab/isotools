// $Id: Hist.cpp,v 1.2 2011/01/08 01:13:10 samn Exp $ 
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

#include "Hist.h"
#include "Vertex.h"
#include "Cluster.h"

using namespace std;

void GetFullBGDistrib(vector<float>& vFloat,KDTreeHist& oTree,int iDims,int* pBestDims,int iBestDims)
{
	int i = 0, j = 0 , iV = 0 , iTotalVs = vFloat.size() / iDims;
	vector<float> vFullData(iTotalVs*iBestDims);
	for(iV=0;iV<iTotalVs;iV++)
	{	for(i=0;i<iBestDims;i++)
			vFullData[j++]=vFloat[iV*iDims+pBestDims[i]];
	}
	oTree.SetData(iBestDims,&vFullData[0],iTotalVs);
}

//this is the continuous multidimensional probability version
void FillDistribs(vector<float>& vFloat,vector<KDTreeHist>& vDistribs,vector<KDTreeHist>& vCompDistribs,int iDistribs,vector<int>& vClustIDs,vector<int>& vCounts,int iDims,A2D<int>& vBestDims,int iBestDims)
{
	vDistribs = vector< KDTreeHist >(iDistribs+1);
	vCompDistribs = vector< KDTreeHist >(iDistribs+1);

	int iTotalVs = vFloat.size() / iDims, iC = 0;

	//full distribution not really used so no need to initialize!!!!!

	for(iC=1;iC<iDistribs;iC++)
	{	int iCompSize = iTotalVs - vCounts[iC];
		vector<float> vClustData(vCounts[iC]*iBestDims), vCompData(iCompSize*iBestDims);
		int i = 0, j = 0 , k = 0, iV = 0;
		for(iV=0;iV<vClustIDs.size();iV++)
		{	if(vClustIDs[iV] == iC)
			{	for(i=0;i<iBestDims;i++)
					vClustData[j++]=vFloat[iV*iDims+vBestDims[iC][i]];
			}
			else 
			{	for(i=0;i<iBestDims;i++)
					vCompData[k++]=vFloat[iV*iDims+vBestDims[iC][i]];
			}
		}
		vDistribs[iC].SetData(iBestDims,&vClustData[0],vCounts[iC]);
		vCompDistribs[iC].SetData(iBestDims,&vCompData[0],iCompSize);
	}
}

static vector< vector<float> > gvprobs;
void InitProbs(int iMaxNumElems)
{
	if(gvprobs.size()>=iMaxNumElems+1)return;
	gvprobs = vector< vector<float> >(iMaxNumElems+1);
	int i,j;
	gvprobs[0] = vector<prob_t>(1);
	gvprobs[0][0] = 0.0;
	for(i=1;i<=iMaxNumElems;i++)
	{
		gvprobs[i] = vector<prob_t>(i+1);
		for(j=0;j<=i;j++)
		{
			gvprobs[i][j] = (prob_t) j / (prob_t) i;
		}
	}
}

float Prob(int iElems,int i)
{
	if(iElems >= gvprobs.size() || i >= gvprobs[iElems].size())
		return (float) i / (float) iElems;
	return gvprobs[iElems][i];
}

ProbInitFree::ProbInitFree(int i)
{
	InitProbs(i);
}

ProbInitFree::~ProbInitFree()
{
	gvprobs.clear();
}
