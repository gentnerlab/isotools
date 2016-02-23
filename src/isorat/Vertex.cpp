// $Id: Vertex.cpp,v 1.4 2011/01/08 01:15:04 samn Exp $ 
// Vertex.cpp: implementation of the CVertex class.
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

#include "Vertex.h"
#include "Log.h"
#include "WCMath.h"
#include "StringUtils.h"
#include <algorithm>
#include <math.h>
#include <errno.h> 
#include <fstream>

/////////////////////////////////////////////////////////////////////////////

void CVertex::AddPnt(float toStore)
{
	m_Vertex.push_back(toStore);
}

/////////////////////////////////////////////////////////////////////////////
// CVertexStack
CVerxStack::CVerxStack()
{ 
	m_NumVerx = 0;
	m_iNumDim = 0;
	m_bNormFloatV = true;
}

void CVerxStack::CalcDimStats()
{
	MY_STACK::iterator iVerx;	
	m_MainMin.clear();
	m_MainMax.clear();
	m_MainRange = VERTEX(m_iNumDim,0.0f);
	int i;
	for (i=0;i<m_iNumDim;i++)
	{	m_MainMin.push_back((float) 2e+20);
		m_MainMax.push_back((float)-2e+20);
	}	
	for (iVerx = m_VerxStack.begin(); iVerx != m_VerxStack.end(); iVerx++)
	{	CVertex* verx = (CVertex*) *iVerx;		
		for (i=0;i<m_iNumDim;i++)			
		{	float val = verx->GetValue(i);
			if ( val < m_MainMin[i] )
				m_MainMin[i] = val;
			if ( val > m_MainMax[i] )
				m_MainMax[i] = val;
		}
	}
	for(i=0;i<m_iNumDim;i++) 
		m_MainRange[i] = m_MainMax[i] - m_MainMin[i];
}

/////////////////////////////////////////////////////////////////////////////
void CVerxStack::CalcMinMax()
{
	MY_STACK::iterator iVerx;	
	m_MainMin.clear();
	m_MainMax.clear();
	int i;
	for (i=0;i<m_iNumDim;i++)
	{
		m_MainMin.push_back((float) 2e+20);
		m_MainMax.push_back((float)-2e+20);
	}	
	CVertex *verx;
	for (iVerx = m_VerxStack.begin(); iVerx != m_VerxStack.end(); iVerx++)
	{
	  verx = (CVertex*) *iVerx;	 
	  for (i=0;i<m_iNumDim;i++)
	  {
	  	if ( verx->GetValue(i) < m_MainMin[i])
	  		m_MainMin[i] = verx->GetValue(i);
	    if ( verx->GetValue(i) > m_MainMax[i])
	    	m_MainMax[i] = verx->GetValue(i);
	  }
	}
}

/////////////////////////////////////////////////////////////////////////////
float CVerxStack::GetMin(int Index)
{ 
  if (Index<m_iNumDim) 
    return m_MainMin[Index];
  else 
    return -10; 
}

/////////////////////////////////////////////////////////////////////////////
float CVerxStack::GetMax(int Index)
{ 
  if (Index<m_iNumDim)
    return m_MainMax[Index];
  else 
    return 10;
}

/////////////////////////////////////////////////////////////////////////////
void CVerxStack::AddVrx(CMyObject *toStore)
{
  m_VerxStack.push_back(toStore);
}

/////////////////////////////////////////////////////////////////////////////
int CVerxStack::GetClust(int NoOfVerx)
{
  if (NoOfVerx < m_NumVerx)
    {
      return ((CVertex*)*(m_VerxStack.begin() + NoOfVerx))->GetClust();
    }
  return -1;
}

void InitProbs(int iMaxNumElems);

/////////////////////////////////////////////////////////////////////////////
void CVerxStack::CalcAfterLoad()
{	
//#ifdef _DEBUG
  CalcMinMax();
  //int i;
  //for(i=0;i<m_iNumDim;i++) printf("%g %g\n",m_MainMin[i],m_MainMax[i]);
  //WriteVec2Log(this->m_MainMin);	
  //WriteVec2Log(this->m_MainMax);  
//#endif
}

bool CVerxStack::ReadAsciiData(char* path) {	
  try {
    int i = 0, iMaxClust = 0;
    char buf[2048];
    vector<string> vstr;
    string delim(" \t");
    // do actual loading here
    // line 1 is header (space-separated column names, column 0 is cluster id)
    // line 2... n - vector record
    //record structure : clusterID feature1 ... featuren
    std::ifstream f;
    f.open(path);
    if(!f.is_open()) return false;      
    f.getline(buf,2048);
    Split(buf,delim,vstr);
    m_iNumDim = vstr.size() - 1; // # of dimensions (first column is cluster ID)
    if(m_iNumDim < 2) {
    	fprintf(stderr,"Must have at least 2 dimensions!\n");
    	return false;    
    }
    vector<float> v(m_iNumDim,0.0);
    while(true) {
      int cid;
      f >> cid;      
      if(f.eof()) break;      
      for(i=0;i<m_iNumDim && !f.eof();i++) f >> v[i];      
      if(i<m_iNumDim && f.eof()) break; // full vector not read in 		
	  if(cid > iMaxClust) iMaxClust = cid;         
      CVertex* NewVerx=new CVertex();// Storing data to container
      NewVerx->SetClust(cid); 
      NewVerx->m_Vertex.resize(m_iNumDim);
      copy(v.begin(),v.end(),NewVerx->m_Vertex.begin());      
      m_VerxStack.push_back(NewVerx);		
      if(f.eof()) break;
    }
    f.close();    
    m_NumVerx = m_VerxStack.size();
    if (m_NumVerx == 0) {//check if spikes were loaded
      	return false;
    } else {
    	InitClust(iMaxClust);
      	CalcAfterLoad();
      	Write2Log("Imported %d spikes",m_NumVerx);
      	return true;
    }
  } catch(...) {
    	Write2Log("Exception in ReadData!");
    	return false;
  }
}

void CVerxStack::InitClust(int iMaxClust)
{
	m_oClusters.m_iNumClusts =iMaxClust; 
}

/////////////////////////////////////////////////////////////////////////////
void CVerxStack::SetEmpty()
{
	  // Remove vectors
  	MY_STACK::iterator Index;
  	for (Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++)
    {
	      CVertex* verx = (CVertex*)*Index;
    	  delete verx;
    }
  	m_VerxStack.clear();
  	m_iNumDim = 0;
  	m_NumVerx = 0;
}
