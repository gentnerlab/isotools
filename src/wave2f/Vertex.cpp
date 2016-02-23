// $Id: Vertex.cpp,v 1.13 2011/07/06 04:13:10 samn Exp $ 
// Vertex.cpp: implementation of the CVertex class.
//
/*  wave2f - This program calculates the features from a text file consisting
    of waveforms from neural extracellular field recordings using tetrodes.
    The features include peak, energy, and optionally, principal components
    of the waveforms.

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
*/

//////////////////////////////////////////////////////////////////////

#include "Vertex.h"
#include "pca2.hpp"
#include <algorithm>
#include <math.h>
#include <errno.h> 
#include <fstream>

/////////////////////////////////////////////////////////////////////////////

int CVertex::FREQ = 10000;
int CVertex::SAMPLES = 64;
int CVertex::RESOLUTION = 16;
int CVertex::CHANNELS = 4;
vector<float> CVertex::VX;

void CVertex::AddPnt(float toStore) // adds a feature to m_Vertex
{
  m_Vertex.push_back(toStore);
}

prob_t CVertex::GetEnergy(int iSamples,int iChannel)
{
  prob_t dEnergy = 0.0;
  float* y = &m_YValues[iSamples*iChannel];
  int i;
  for(i=0;i<iSamples;i++)
    dEnergy += y[i]*y[i];
  dEnergy = sqrt(dEnergy);
  if(iSamples) 
    dEnergy /= iSamples;
  return dEnergy;
}

/////////////////////////////////////////////////////////////////////////////
void CVertex::SetYValues(vector<short>& vBuffer, vector<float>& x, vector<float>& u, int channels, int samples)
{
  int iSz = vBuffer.size(),i;

  //ensure proper size for buffers
  m_YValues.resize(iSz);
  m_d2YValues.resize(iSz);
	
  // Convert to real value
  for (int m_l=0;m_l<iSz;m_l++)
    m_YValues[m_l] = (float) vBuffer[m_l];
  //m_YValues[m_l] = (float) (10 * vBuffer[m_l])/SHRT_MAX;

  //do interpolation
  
  for (i=0; i<channels; i++)
    Spline(&x[0],&m_YValues[samples*i],samples,&m_d2YValues[samples*i],&u[0]);
}

/////////////////////////////////////////////////////////////////////////////
void CVertex::SetYValues(vector<float>& mNew, vector<float>& x,vector<float>& u, int channels, int samples)
{ int i;
  m_YValues.resize(mNew.size());
  copy(mNew.begin(),mNew.end(),m_YValues.begin());
  m_d2YValues.resize(channels*samples);
  for (i=0; i<channels; i++)
    Spline(&x[0],&m_YValues[samples*i],samples,&m_d2YValues[samples*i],&u[0]);
}

/////////////////////////////////////////////////////////////////////////////
float CVertex::GetYValues(int Index)
{ 
  if ( Index>=0 && Index<m_YValues.size()) 
    return m_YValues[Index]; 
  else 
    return 0; 
}

void CVertex::ExportFeatures(FILE* fp)
{
  fprintf(fp,"%d\t", GetClust() ); // cluster ID
  VERTEX::iterator index;
  for (index = m_Vertex.begin(); index != m_Vertex.end(); index++) // m_Vertex stores features
    {
      fprintf(fp,"%f\t",*index); // feature value
    }
  fprintf(fp,"\n");
}

/////////////////////////////////////////////////////////////////////////////
// CVertexStack
CVerxStack::CVerxStack(bool bWCF, int iNumPCTs, int iNumPCCs)
{ 
  m_NumVerx = 0;
  m_iNumDim = 0;
  m_bNormFloatV = true;

  RESOLUTION = 16;
  NUM_CHANNELS = 4; 
  NUM_SAMPLES = 64;
  SAMPLE_FREQ = 10000; 
  m_x.resize(NUM_SAMPLES);
  CVertex::VX.resize(NUM_SAMPLES);
  for (int i=0;i<NUM_SAMPLES;i++)
    CVertex::VX[i] = m_x[i] = (float) i/SAMPLE_FREQ;

  m_bWCF = bWCF;
  m_iNumPCTs = iNumPCTs;
  m_iNumPCCs = iNumPCCs;
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
    { CVertex* verx = (CVertex*) *iVerx;		
      for (i=0;i<m_iNumDim;i++)			
	{ float val = verx->GetValue(i);
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
  MY_STACK::iterator Index;
  for (Index=m_VerxStack.begin(); Index!=m_VerxStack.end(); Index++)
  {
    CVertex* NewVerx = (CVertex*) *Index;
    CalcOneSpike(NewVerx);
  }

//#ifdef _DEBUG
  CalcMinMax();
  //int i;
  //for(i=0;i<m_iNumDim;i++) printf("%g %g\n",m_MainMin[i],m_MainMax[i]);
  //WriteVec2Log(this->m_MainMin);	
  //WriteVec2Log(this->m_MainMax);  
//#endif

  if(m_iNumPCTs > 0 || m_iNumPCCs > 0) CalcPCA();
}

void CVerxStack::CalcPCA()
{
  if(m_iNumPCTs > 0) CalcPCAT();
  if(m_iNumPCCs > 0) CalcPCAC();
}

void CVerxStack::CalcPCAT()
{
  if(m_iNumPCTs < 1) return;
  const int iParams = m_iNumPCTs;
  //vStore = storage for principle components. later will move it to the vertices.
  //only used to keep a specific ordering on the data
  A2D<float> vStore(m_VerxStack.size(),iParams*NUM_CHANNELS);
  int iChannel = 0, iEC = m_bWCF ? 16 : 4; // iEC is energy param index
  for(iChannel=0;iChannel<NUM_CHANNELS;iChannel++) {	
    int iV = 0, iPCASz = NUM_SAMPLES;
    MY_STACK::iterator Index;  
    A2D<PCA_T> vChanData(m_VerxStack.size(),iPCASz);
    PCA oPCA(iPCASz);
    for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++) { //copy vertex data to pca format
      PCA_T* vData = vChanData[iV];
      CVertex* verx = (CVertex*) *Index;
      prob_t dEnergyNorm = verx->GetValue(iEC+iChannel) * NUM_SAMPLES;
      int i = 0;
      if(dEnergyNorm <= 0.0) dEnergyNorm = 1.0;
      for(i=0;i<NUM_SAMPLES;i++) vData[i] = verx->GetYValues((iChannel*NUM_SAMPLES)+i) / dEnergyNorm;
    }
    for(iV=0;iV<vChanData.Rows();iV++)  oPCA.putData(vChanData[iV]);
    oPCA.dataEnd();  oPCA.calcPCA();    
    vector<PCA_T> vTrans(iParams); //store projections of row-points on first 3 prin. comps.
    for(iV=0;iV<m_VerxStack.size();iV++) {
      oPCA.transform(vChanData[iV],iPCASz,&vTrans[0],iParams);
      int iPC = 0;
      for(iPC = 0; iPC < iParams; iPC++) vStore[iV][iChannel*iParams+iPC] = vTrans[iPC];
    }    
  }  
  int iV=0; //store the projection onto the principle components in vertices
  MY_STACK::iterator Index;
  for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++) {
    CVertex* verx = (CVertex*) *Index;  int iP;
    for(iP=0;iP<iParams;iP++) for(iChannel=0;iChannel<NUM_CHANNELS;iChannel++)verx->AddPnt(vStore[iV][iParams*iChannel+iP]);
  }		
}

void CVerxStack::CalcPCAC()
{
  if(m_iNumPCCs < 1) return;
  const int iParams = m_iNumPCTs;
  //vStore = storage for principle components. later will move it to the vertices.
  //only used to keep a specific ordering on the data
  A2D<float> vStore(m_VerxStack.size(),iParams*NUM_CHANNELS);
  int iV = 0, iPCASz = NUM_SAMPLES * NUM_CHANNELS;
  MY_STACK::iterator Index;  
  A2D<PCA_T> vChanData(m_VerxStack.size(),iPCASz);
  PCA oPCA(iPCASz);
  int iSamps = iPCASz;
  vector<double> vtmp(iSamps), vEnergyNorm(m_VerxStack.size());
  for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++) {  //copy vertex data to pca format
    CVertex* verx = (CVertex*) *Index;
    int i = 0;
    for(i=0;i<iSamps;i++) vtmp[i]=verx->GetYValues(i);
    double dEnergyNorm =  Energy(vtmp) * iSamps;
    vEnergyNorm[iV] = dEnergyNorm;
    for(i=0;i<iSamps;i++) vtmp[i] /= dEnergyNorm;
    oPCA.putData(vtmp);
  }
  oPCA.dataEnd(); oPCA.calcPCA(); iV=0;
  vector<double> vout(m_iNumPCCs,0.0);
  for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++,iV++) {
    CVertex* verx = (CVertex*) *Index;
    int iP , i;
    for(i=0;i<iSamps;i++) vtmp[i]=verx->GetYValues(i) / vEnergyNorm[iV];
    fill(vout.begin(),vout.end(),0.0f);
    oPCA.transform(&vtmp[0],iSamps,&vout[0],m_iNumPCCs);
    for(iP=0;iP<m_iNumPCCs;iP++)//store projections of spikes on PCs
      verx->AddPnt(vout[iP]);
  }
}

bool CVerxStack::ImportSpikes(char* path,bool bAscii)
{
  SetEmpty();
  CVertex::CHANNELS = NUM_CHANNELS = 4; 
  CVertex::SAMPLES = NUM_SAMPLES = 32;
  if(bAscii) {
    if(!ImportSpikesAscii(path)) return false;
  } else if(!ImportSpikesBin(path)) return false;
  if (m_NumVerx == 0) return false; //check if spikes were loaded 
  CalcAfterLoad(); // calc waveform features
  return true;
}

bool CVerxStack::ImportSpikesBin(char* path)
{ try{
  // file format:
  // int 1 - # of channels (i.e., 4 for a tetrode)
  // int 2 - # of samples per spike on single wire
  // int 3 - sampling frequency
  // rest of file is waveform records
  //record structure : clusterID(int) numsamples*numwires(double)
  //  samples arranged in order of wires on tetrode (or septode,etc.)
  FILE* fp;
  if(!(fp = fopen(path,"rb"))) return false;  

  vector<int> viBuf(3);
  if( 3 != fread(&viBuf[0], sizeof(int), 3, fp) ) return false;

  NUM_CHANNELS = CVertex::CHANNELS = viBuf[0]; 
  NUM_SAMPLES = CVertex::SAMPLES = viBuf[1];
  int samps = NUM_CHANNELS*NUM_SAMPLES;
  if(!samps) return false;

  m_x.resize(NUM_SAMPLES);
  CVertex::VX.resize(NUM_SAMPLES);
  int i;
  for (i=0;i<NUM_SAMPLES;i++) CVertex::VX[i] = m_x[i] = (float) i/SAMPLE_FREQ;

  SAMPLE_FREQ = CVertex::FREQ = viBuf[2];
  if(SAMPLE_FREQ <= 0.0) return false;

  vector<double> vd(samps,0.0); 
  vector<float> m_u(NUM_SAMPLES) , v(samps,0.0);
  while(true) { int cid;
    if( 1 != fread(&cid,sizeof(int),1,fp) ) break; // cluster ID
    if( samps != fread(&vd[0],sizeof(double),samps,fp) ) break; // read in waveforms
    for(i=0;i<samps;i++) v[i] = (float) vd[i]; // copy to floats
    CVertex* NewVerx=new CVertex(); // Storing data to container
    NewVerx->SetClust(cid);    
    NewVerx->SetYValues(v, m_x, m_u, NUM_CHANNELS, NUM_SAMPLES);
    m_VerxStack.push_back(NewVerx);
    m_NumVerx++;
    if(feof(fp)) break;
  }
  fclose(fp);
  return true;
  } catch(...) {
    return false;
  }
}

bool CVerxStack::ImportSpikesAscii(char* path)
{ try{
    int i = 0;
    // do actual loading here
    // line 1 - # of tetrodes
    // line 2 - # of samples per spike on single wire
    // line 3 - sampling frequency
    // line 4... n - waveform record
    //record structure : clusterID numsamples*numtetrodes
    //first waveform for tetrode 1 , then tetrode 2 , then 3, then 4
    std::ifstream f;
    f.open(path);
    if(!f.is_open())
      return false;
    f >> CVertex::CHANNELS;
    NUM_CHANNELS = CVertex::CHANNELS;
    f >> CVertex::SAMPLES;
    NUM_SAMPLES = CVertex::SAMPLES;
    int samps = NUM_CHANNELS*NUM_SAMPLES;
    if(!samps) return false;    

    m_x.resize(NUM_SAMPLES);
    CVertex::VX.resize(NUM_SAMPLES);
    for (i=0;i<NUM_SAMPLES;i++) CVertex::VX[i] = m_x[i] = (float) i/SAMPLE_FREQ;

    vector<float> v(samps,0.0) , m_u(NUM_SAMPLES);
    f >> CVertex::FREQ;
    SAMPLE_FREQ = CVertex::FREQ;
    if(SAMPLE_FREQ <= 0.0) return false;
    while(true) {
      int cid;
      f >> cid; // cluster ID
      if(f.eof()) break;
      for(i=0;i<samps && !f.eof();i++) f >> v[i]; // read in features
      if(i<samps && f.eof()) break; // full waveform not read in
      CVertex* NewVerx=new CVertex(); // Storing data to container
      NewVerx->SetClust(cid);									
      m_NumVerx++;	
      NewVerx->SetYValues(v, m_x, m_u, NUM_CHANNELS, NUM_SAMPLES);	
      m_VerxStack.push_back(NewVerx);		
      if(f.eof()) break;
    }
    f.close(); 
    return true;
  }
  catch(...){
    return false;
  }
}

bool CVerxStack::ReadAsciiData(char* path) {	
  try {
    int i = 0, iMaxClust = 0;
    // do actual loading here
    // line 1 - # of dimensions 
    // line 2... n - vector record
    //record structure : clusterID feature1 ... featuren
    std::ifstream f;
    f.open(path);
    if(!f.is_open()) return false;      
    f >> m_iNumDim; // # of dimensions  
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
      	CalcAfterLoad();
      	return true;
    }
  } catch(...) {
    	return false;
  }
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

// Written by Andre Fenton
void CVertex::Spline(float *x,float *y,int n,float *d2y,float *u)
{
	register int	i;
	register float	p, sig;

	d2y[0] = u[0] = 0.0;

	for (i = 1; i < n - 1; i++) {
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig * d2y[i - 1] + 2.0;
		d2y[i] = (sig - 1.0) / p;
		u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
		u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
	}

	d2y[n - 1] = 0.0;

	for (i = n - 1; i--; )
		d2y[i] = d2y[i] * d2y[i + 1] + u[i];

	return;
}

void CVertex::Spline(float *x,int n,float *u)
{
	register int	i;
	register float	p, sig;
	float *y,*d2y;
	y = &m_YValues[0];
	d2y = &m_d2YValues[0];

	d2y[0] = u[0] = 0.0;

	for (i = 1; i < n - 1; i++) {
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig * d2y[i - 1] + 2.0;
		d2y[i] = (sig - 1.0) / p;
		u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
		u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
	}

	d2y[n - 1] = 0.0;

	for (i = n - 1; i--; )
		d2y[i] = d2y[i] * d2y[i + 1] + u[i];

	return;
}


/////////////////////////////////////////////////////////////////////////////
// Written by Andre Fenton
float CVertex::Splope(float *xa,float *ya,float *d2y,float x,int lo_pt,int hi_pt)
{
	register float a, b, c, d, e, f, g, h;

	a = xa[lo_pt];
	b = xa[hi_pt];
	c = ya[lo_pt];
	d = ya[hi_pt];
	e = d2y[lo_pt];
	f = d2y[hi_pt];
	g = b - a;

	h = (x * (e * b - f * a) + d - c) / g;
	h += (x * x * (f - e) + f * a * a - e * b * b) / (2.0 * g);
	h += g * (e - f) / 6.0;

	return(h);
}

float CVertex::Splope(float *xa,float x,int lo_pt,int hi_pt)
{
	float *ya,*d2y;
	ya = &m_YValues[0];
	d2y = &m_d2YValues[0];
	
	register float a, b, c, d, e, f, g, h;

	a = xa[lo_pt];
	b = xa[hi_pt];
	c = ya[lo_pt];
	d = ya[hi_pt];
	e = d2y[lo_pt];
	f = d2y[hi_pt];
	g = b - a;

	h = (x * (e * b - f * a) + d - c) / g;
	h += (x * x * (f - e) + f * a * a - e * b * b) / (2.0 * g);
	h += g * (e - f) / 6.0;

	return(h);
}

float CVertex::Splope(float *xa,float x,int lo_pt,int hi_pt,int channel,int Samples)
{
	float *ya,*d2y;
	ya = &m_YValues[channel*Samples];
	d2y = &m_d2YValues[channel*Samples];
	
	register float a, b, c, d, e, f, g, h;

	a = xa[lo_pt];
	b = xa[hi_pt];
	c = ya[lo_pt];
	d = ya[hi_pt];
	e = d2y[lo_pt];
	f = d2y[hi_pt];
	g = b - a;

	h = (x * (e * b - f * a) + d - c) / g;
	h += (x * x * (f - e) + f * a * a - e * b * b) / (2.0 * g);
	h += g * (e - f) / 6.0;

	return(h);
}

/////////////////////////////////////////////////////////////////////////////
// Written by Andre Fenton
float CVertex::Splint(float *xa,float *ya,float *d2y,float x,int lo_pt,int hi_pt)
{
	register float	h, b, a;
	float	y;

	h = xa[hi_pt] - xa[lo_pt];
	a = (xa[hi_pt] - x) / h;
	b = (x - xa[lo_pt]) / h;

	y = a * ya[lo_pt] + b * ya[hi_pt] + ((a * a * a - a) * d2y[lo_pt] + (b * b * b - b) * d2y[hi_pt]) * (h * h) / 6.0;

	return(y);
}

float CVertex::Splint(float *xa,float x,int lo_pt,int hi_pt)
{
	float *ya,*d2y;
	ya = &m_YValues[0];
	d2y = &m_d2YValues[0];
		
	register float	h, b, a;
	float	y;

	h = xa[hi_pt] - xa[lo_pt];
	a = (xa[hi_pt] - x) / h;
	b = (x - xa[lo_pt]) / h;

	y = a * ya[lo_pt] + b * ya[hi_pt] + ((a * a * a - a) * d2y[lo_pt] + (b * b * b - b) * d2y[hi_pt]) * (h * h) / 6.0;

	return(y);
}

float CVertex::Splint(float *xa,float x,int lo_pt,int hi_pt,int channel,int Samples)
{
	float *ya,*d2y;
	ya = &m_YValues[channel*Samples];
	d2y = &m_d2YValues[channel*Samples];
		
	register float	h, b, a;
	float	y;

	h = xa[hi_pt] - xa[lo_pt];
	a = (xa[hi_pt] - x) / h;
	b = (x - xa[lo_pt]) / h;

	y = a * ya[lo_pt] + b * ya[hi_pt] + ((a * a * a - a) * d2y[lo_pt] + (b * b * b - b) * d2y[hi_pt]) * (h * h) / 6.0;

	return(y);
}

void CVerxStack::SplineUpSample(CVertex* verx,float* vout,int iChannel,int iINC)
{
  int iK=0 , idx = 0, iL = 0;
  for (iK=0;iK<NUM_SAMPLES-1;++iK)
    for (iL=0; iL<RESOLUTION; iL+=iINC)
      vout[idx++]=verx->Splint(&m_x[0], m_x[iK]+(float)iL/(SAMPLE_FREQ*RESOLUTION), iK+0, iK+1, iChannel, NUM_SAMPLES);
  vout[idx++]=verx->GetYValues(iChannel*NUM_SAMPLES+NUM_SAMPLES-1);
}

int CVerxStack::MaxIndex(CVertex* verx,int iChannel,int iINC)
{
  int m_iK;
  int lo_ptMax;
  float m_max,mTMax;	
  m_max=(float)-2e+20;						
  int iQ = iChannel;							
  m_iK=0;
  float dYt1 = verx->Splope(&m_x[0], m_x[m_iK], m_iK+0, m_iK+1, iQ, NUM_SAMPLES);
  for (m_iK=0;m_iK<(NUM_SAMPLES-2);++m_iK)
    {
      float dYt2 = verx->Splope(&m_x[0], m_x[m_iK+1], m_iK+1, m_iK+2, iQ, NUM_SAMPLES);
      if (dYt1*dYt2<=0)
	{
	  for (int iL=0; iL<RESOLUTION; iL+=iINC)
	    {
	      float Value = verx->Splint(&m_x[0], m_x[m_iK]+(float)iL/(SAMPLE_FREQ*RESOLUTION), m_iK+0, m_iK+1, iQ, NUM_SAMPLES);
	      if (m_max<Value)
		{
		  m_max=Value;
		  mTMax=m_x[m_iK]+(float)iL/(SAMPLE_FREQ*RESOLUTION);
		  lo_ptMax=m_iK*RESOLUTION/iINC + iL/iINC;
		}
	    } // 
	} // If is extreme, change resolution
      dYt1 = dYt2;
    } // For each 64 (32,...) samples
  return lo_ptMax;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CVerxStack::CalcOneSpike(CVertex *verx) // calculates features for waveforms associated with 1 spike
{
  int m_iK, lo_ptMin[4], lo_ptMax[4], lo_ptMinDef,lo_ptMaxDef, iQ;
  float m_min[4],m_max[4],mTMin[4],mTMax[4],m_minDef,m_maxDef;	
  for (iQ=0;iQ<NUM_CHANNELS;iQ++) {
    m_min[iQ]=(float) 2e+20;
    m_max[iQ]=(float)-2e+20;
  }						
  float dYt1,dYt2;							
  for (iQ=0;iQ<NUM_CHANNELS;iQ++) {
    m_iK=0;
    dYt1 = verx->Splope(&m_x[0], m_x[m_iK], m_iK+0, m_iK+1, iQ, NUM_SAMPLES);
    for (m_iK=0;m_iK<(NUM_SAMPLES-2);++m_iK) {
      dYt2 = verx->Splope(&m_x[0], m_x[m_iK+1], m_iK+1, m_iK+2, iQ, NUM_SAMPLES);
      if (dYt1*dYt2<=0) {
	for (int iL=0; iL<RESOLUTION; iL++) {
	  float Value = verx->Splint(&m_x[0], m_x[m_iK]+(float)iL/(SAMPLE_FREQ*RESOLUTION), m_iK+0, m_iK+1, iQ, NUM_SAMPLES);
	  if (Value >= m_max[iQ] || Value <= m_min[iQ]) {
	    if (m_min[iQ]>Value) {
	      m_min[iQ]=Value;
	      mTMin[iQ]=m_x[m_iK]+(float)iL/(SAMPLE_FREQ*RESOLUTION);
	      lo_ptMin[iQ]=m_iK;
	    }
	    if (m_max[iQ]<Value) {
	      m_max[iQ]=Value;
	      mTMax[iQ]=m_x[m_iK]+(float)iL/(SAMPLE_FREQ*RESOLUTION);
	      lo_ptMax[iQ]=m_iK;
	    }
	  }
	} // 
      } // If is extreme, change resolution
      dYt1 = dYt2;
    } // For each 64 (32,...) samples
  } //  Repeat for each channel  
  
  float MainTMin = mTMin[0];
  float MainTMax = mTMax[0];

  m_minDef = m_min[0];
  m_maxDef = m_max[0];
  lo_ptMinDef = lo_ptMin[0];
  lo_ptMaxDef = lo_ptMax[0];

  for (int iN=1;iN<NUM_CHANNELS;iN++) {
    if (m_min[iN]<m_minDef) { 
      m_minDef = m_min[iN];
      lo_ptMinDef = lo_ptMin[iN];
      MainTMin = mTMin[iN];
    }
    if (m_max[iN]>m_maxDef) { 
      m_maxDef = m_max[iN];
      lo_ptMaxDef = lo_ptMax[iN];
      MainTMax = mTMax[iN];
    }
  }  
  // store features
  float minToStore,maxToStore;  int iP=0;
  for (iP=0;iP<NUM_CHANNELS;iP++) verx->AddPnt(m_max[iP]); // Tx_peak

  if(m_bWCF) { // WClust-style features
    for (iP=0;iP<NUM_CHANNELS;iP++) { 	
      maxToStore = verx->Splint(&m_x[0], MainTMax, lo_ptMaxDef, lo_ptMaxDef+1, iP, NUM_SAMPLES);
      verx->AddPnt(maxToStore);	// Tx-V(peak)
    }  
    for (iP=0;iP<NUM_CHANNELS;iP++) verx->AddPnt(m_min[iP]); // Tx-valley  
    for (iP=0;iP<NUM_CHANNELS;iP++) {
      minToStore = verx->Splint(&m_x[0], MainTMin, lo_ptMinDef, lo_ptMinDef+1, iP, NUM_SAMPLES);
      verx->AddPnt(minToStore);	// Tx-V(valley)
    }
  }    
  for(iP=0;iP<NUM_CHANNELS;iP++) { //store energy
    prob_t dEnergy = verx->GetEnergy(NUM_SAMPLES,iP);
    verx->AddPnt(dEnergy);
  }
}

/////////////////////this header is for when exporting features//////////////
void CVerxStack::WriteHeader(FILE* fp)
{
  vector<string> vs; unsigned int i; int j; char tmp[256];
  vs.push_back(string("Peak"));
  if(m_bWCF) {
    vs.push_back(string("V(peak)"));
    vs.push_back(string("Valley"));
    vs.push_back(string("V(valley)"));
  }
  vs.push_back(string("Energy"));
  for(i=1;i<=m_iNumPCTs;i++) {sprintf(tmp,"PC%d",i); vs.push_back(string(tmp));} // PCA features
  fprintf(fp,"Cluster ");
  for(i=0;i<vs.size();i++)
    for(j=1;j<=NUM_CHANNELS;j++) 
      fprintf(fp,"T%d-%s ",j,vs[i].c_str());
  for(i=1;i<=m_iNumPCCs;i++) fprintf(fp,"PCC%d ",i);
  fprintf(fp,"\n");
}

bool CVerxStack::ExportFeatures(char* path)
{
  FILE* fp = fopen(path,"w");
  if(!fp) return false;
  WriteHeader(fp);
  MY_STACK::iterator iVerx;
  for (iVerx = m_VerxStack.begin(); iVerx != m_VerxStack.end(); iVerx++) { 
    CVertex* verx = (CVertex*) *iVerx;
    verx->ExportFeatures(fp);
  }
  fclose(fp);
  return true;
}

bool CVerxStack::ExportUpsampSpikes(char* path)
{
  FILE* fp = fopen(path,"w");
  if(!fp) return false;
  MY_STACK::iterator Index;
  int iResRedFctr = 4; int i; unsigned int j;
  vector<float> vUp(NUM_SAMPLES*RESOLUTION/iResRedFctr+1,0.0f);//vector storing upsampled spike
  fprintf(fp,"Cluster ");
  for(i=0;i<NUM_CHANNELS;i++)//write header for NQS (or other, can get rid of before load into matlab)
    for(j=0;j<vUp.size();j++)
      fprintf(fp,"C%dS%d ",i,j);
  fprintf(fp,"\n");
  for(Index=m_VerxStack.begin();Index!=m_VerxStack.end();Index++) {	
    CVertex* verx = (CVertex*) *Index;
    fprintf(fp,"%d ",verx->GetClust());//first column cluster
    for(i=0;i<NUM_CHANNELS;i++) { //rest of columns spike waveforms
      SplineUpSample(verx,&vUp[0],i,iResRedFctr);
      for(j=0;j<vUp.size();j++) fprintf(fp,"%g ",vUp[j]);
    }
    fprintf(fp,"\n");
  }		
  fclose(fp);
  return true;
}
