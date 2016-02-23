// $Id: main.cpp,v 1.9 2011/02/01 00:51:09 samn Exp $ 

/*  wave2f - This program calculates the features from a text file consisting
    of waveforms from neural extracellular field recordings using tetrodes,
    and then saves the features to an output file. The features include peak,
    energy, and optionally, principal components of the waveforms.

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

#include "Vertex.h"
#include <stdio.h>
#include <string.h>

void prhelp() {

  fprintf(stderr,"\nwave2f: Copyright (C) 2003-2011 Sam Neymotin & BioSignal Group.\n\n");

  fprintf(stderr,"wave2f - This program calculates the features from a text file consisting\n");
  fprintf(stderr,"of waveforms from neural extracellular field recordings using tetrodes,\n");
  fprintf(stderr,"and then saves the features to an output file. The features include peak,\n");
  fprintf(stderr,"energy, and optionally, principal components of the waveforms.\n\n");

  fprintf(stderr,"Command-line usage: wave2f input output [-u upsampfile][-bin][-wc][-pca n][-pcac n]\n");
  fprintf(stderr," input is the input file path (format described below)\n");
  fprintf(stderr," output is the output file path (format described below)\n");
  fprintf(stderr," -u upsample : upsamples the waveforms using splines and saves to text file upsampfile\n");
  fprintf(stderr," -bin : specifies that input file is in binary format\n");
  fprintf(stderr," -wc : use same features as used in WClust spike sorting program (default off)\n");
  fprintf(stderr," -pca n : calculate/save the first n principal components on the individual wires of the tetrode\n");
  fprintf(stderr," -pcac n : calculate/save the first n principal components of waveform signal\n");
  fprintf(stderr,"           concatenated from all tetrode wires\n\n");

  fprintf(stderr,"To run the program and obtain the waveform features you will need a data\n");
  fprintf(stderr,"input file in either ASCII/text or binary format (both described next).\n\n"); 

  fprintf(stderr,"Input format (ASCII / text ):\n");
  fprintf(stderr," line 1 : number of channels\n");
  fprintf(stderr," line 2 : number of samples on each channel\n");
  fprintf(stderr," line 3 : sampling frequency in Hz\n");
  fprintf(stderr," line 4 through end : waveform record\n");
  fprintf(stderr,"  waveform record is: cluster id followed\n");
  fprintf(stderr,"   by number samples * number channels waveform signal values\n");
  fprintf(stderr,"   each waveform record is on a separate line.\n\n");

  fprintf(stderr,"Input format (BINARY):\n");
  fprintf(stderr," integer : # of channels (4 for tetrode)\n");
  fprintf(stderr," integer : # of samples per spike on single wire\n");
  fprintf(stderr," integer : sampling frequency in Hz\n");
  fprintf(stderr," rest of file is waveform records in binary format\n");
  fprintf(stderr,"  waveform record: clusterID as integer, followed by numbersamples*numberchannels as doubles\n\n");

  fprintf(stderr,"Output format (ASCII / text ):\n");
  fprintf(stderr," line 1 : feature names, separated by spaces\n");
  fprintf(stderr," line 2 through end : feature records\n");
  fprintf(stderr,"  feature record is: cluster_identifier feature_1 feature_2 ... feature_N\n");
  fprintf(stderr,"  values in feature record are separated by tabs (\t), each feature record is on a separate line\n\n");

  fprintf(stderr,"Once the feature file is generated, it can be used in the isoi and isorat programs\n");
  fprintf(stderr,"to calculate the cluster quality values.\n\n");

}

int main (int argc, char** argv) {
  if(argc < 3) {prhelp();  return 0;}
  bool bAscii = true, bUpSample = false, bWCF = false;
  int iArg = 3, iNumPCTs = 0, iNumPCCs = 0; char* fUpName = 0;
  while(iArg<argc) {
    if(!strcmp(argv[iArg],"-u")) {
      if(iArg+1>=argc) {prhelp(); return 0;}
      iArg++;
      fUpName = argv[iArg];
      bUpSample = true;
      iArg++;
    }
    else if(!strcmp(argv[iArg],"-bin")){bAscii=false; iArg++;}
    else if(!strcmp(argv[iArg],"-wc")){bWCF=true; iArg++;}
    else if(!strcmp(argv[iArg],"-pca")) {
      if(iArg+1>=argc) {prhelp(); return 0;}
      iArg++;
      iNumPCTs = atoi(argv[iArg]);
      if(iNumPCTs < 0) {prhelp(); return 0;}
      iArg++;
    }
    else if(!strcmp(argv[iArg],"-pcac")) {
      if(iArg+1>=argc) {prhelp(); return 0;}
      iArg++;
      iNumPCCs = atoi(argv[iArg]);
      if(iNumPCCs < 0) {prhelp(); return 0;}
      iArg++;
    }
    else {prhelp(); return 0;}
  }
  CVerxStack cs(bWCF,iNumPCTs,iNumPCCs);
  if(cs.ImportSpikes(argv[1],bAscii));
  // printf("read %d %d-dim vectors from %s\n",cs.m_NumVerx,cs.m_iNumDim,argv[1]);
  else {
  	fprintf(stderr,"wave2f ERR0: couldn't read vectors from %s!\n",argv[1]);
  	return 0;
  }  
  int iRet = 1;
  if(cs.ExportFeatures(argv[2])) {
    printf("exported features to %s\n",argv[2]);
  } else {
    fprintf(stderr,"wave2f ERR1: couldn't open output file %s!\n",argv[2]);
    iRet = 0;
  }

  if(bUpSample) {
    if(cs.ExportUpsampSpikes(fUpName))
      printf("exported upsampled spikes to %s\n",fUpName);
    else {
      iRet = 0;
      fprintf(stderr,"wave2f ERR2: couldn't export upsampled spikes to %s\n",fUpName);
    }
  }
  return iRet;
}
