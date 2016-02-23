
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

This program was tested/developed on LINUX systems, but can be compiled to
run on Microsoft Windows or Mac OS. To compile, you will need a standard C++
compiler and the make utility. compile the C++ files from the command line with
the make command. That will produce an architecture-dependent program called isoi.
The intermediate .o files can be removed.

Command-line usage: isoi input output [-wc][-time]

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
 Each line consists of: IsoI_BG (cluster separation from background),IsoI_NN (cluster separation
 from nearest neighboring cluster), ClosestCluster (nearest neighboring cluster identifier).

 The -wc flag specifies that input features are the same as those used in the WClust spike-sorting
 program. When this flag is provided, 2D slices consisting of highly correlated dimensions, such as
 T1-Peak and T1-V(peak) are discarded. If you are not using the -wc flag, it's up to you to ensure
 the dimensions provided to isoi are not overly correlated.

 The -time flag specifies to perform timing test, if set, will add time column for each cluster.

To run the program and obtain the IsoI measures you will need a data
input file in ASCII format that contains multidimensional vectors with
cluster identifiers for each data vector. The first line of input is the
# of dimensions in each vector. The following lines are the records
for each vector. Each record consists of: clusterid value1 ... valueN
followed by a newline. clusterid is an integer, while the values are 
floating point numbers in decimal notation. The isoi program writes one
line of output for each cluster to a text file. 

For more information contact Sam Neymotin ( samn at neurosim dot downstate dot edu )
 or Andre Fenton ( afenton at nyu dot edu ).

References:
 The methods used in this program are described in an article in press
  at The Journal of Neuroscience,
  Measuring the quality of neuronal identification in ensemble recordings
  by Neymotin SA, Lytton WW, Olypher AO, Fenton AA (2011).

