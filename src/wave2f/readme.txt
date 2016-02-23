
wave2f - This program calculates the features from a text file consisting
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

This program was tested/developed on LINUX systems, but can be compiled to
run on Microsoft Windows or Mac OS. To compile, you will need a standard C++
compiler and the make utility. Compile the C++ files from the command line
with the make command. That will produce an architecture-dependent program
called wave2f (short-hand for waveforms to features). The intermediate
.o files can be removed.

Command-line usage: wave2f input output [-u upsampfile][-bin][-wc][-pca n][-pcac n]
 input is the input file path (format described below)
 output is the output file path (format described below)
 -u upsample : upsamples the waveforms using splines and saves to text file upsampfile
 -bin : specifies that input file is in binary format
 -wc : use same features as used in WClust spike sorting program (default off)
 -pca n : calculate/save the first n principal components on the individual wires of the tetrode
 -pcac n : calculate/save the first n principal components of waveform signal
           concatenated from all tetrode wires

To run the program and obtain the waveform features you will need a data
input file in either text or binary format (described next). 

Input format (ASCII / text ):
 line 1 : number of channels
 line 2 : number of samples on each channel
 line 3 : sampling frequency in Hz
 line 4 through end : waveform record
  waveform record is: cluster id followed
   by number samples * number channels waveform signal values
   each waveform record is on a separate line.

Input format (BINARY):
 integer : # of channels (4 for tetrode)
 integer : # of samples per spike on single wire
 integer : sampling frequency in Hz
 rest of file is waveform records in binary format
 waveform record: clusterID as integer, followed by numbersamples*numberchannels as doubles

Output format (ASCII / text ):
 line 1 : feature names, separated by spaces
 line 2 through end : feature records
  feature record is: cluster_identifier feature_1 feature_2 ... feature_N
  values in feature record are separated by tabs (\t), each feature record is on a separate line

Once the feature file is generated, it can be used in the isoi and isorat programs
to calculate the cluster quality values.

For more information contact Sam Neymotin ( samn at neurosim dot downstate dot edu )
 or Andre Fenton ( afenton at nyu dot edu ).

References:
 Some of the methods used in this program are described in an article in press
  at The Journal of Neuroscience,
  Measuring the quality of neuronal identification in ensemble recordings
  by Neymotin SA, Lytton WW, Olypher AO, Fenton AA (2011).

