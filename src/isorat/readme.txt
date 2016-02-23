
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

This program was tested/developed on LINUX systems, but can be compiled to
run on Microsoft Windows or Mac OS. To compile, you will need a standard C++
compiler and the make utility.

compile the C++ files from the command line with the make command.
That will produce an architecture-dependent program called isorat 
(short-hand for Isolation Distance and L-Ratio).

Command line usage: isorat input output [-time]

 The input file is in ASCII/text format.
 The first line of input is a header for the names of the dimensions on the following lines.
 The dimension names are separated by white-space.
 Note that the first dimension must be the cluster identifier.
 After the header, the remaining lines of the input file are feature vector records.
 Each feature vector record consists of: clusterid as an integer followed by decimal values
 of the feature vectors, separated by spaces or tabs. Each record is on a separate line.
 NOTE that a clusterid of 0 indicates the noise distribution. Clusters that correspond
 to isolated single units (non-noise) should have a clusterid greater or equal to 1.

 The output is saved as a text file with one line for each input cluster.
 Each line consists of: IsoD (cluster separation from background)
 and LRatio (cluster separation from background).

 The -time flag specifies to perform timing test, if set, will add time column for each cluster.

For more information contact Sam Neymotin ( samn at neurosim dot downstate dot edu )
 or Andre Fenton ( afenton at nyu dot edu ).

References:
 The methods used in this program are described in an article under review
  at Journal of Neuroscience,
  Measuring the quality of neuronal identification in ensemble recordings
  by Neymotin SA, Lytton WW, Olypher AO, Fenton AA (2011).

  and in, Quantitative measures of cluster quality for use in extracellular recordings
  Schmitzer-Torbert N, Jackson J, Henze D, Harris K, Redish AD
  Neuroscience 131(1):1-11, 2005.
