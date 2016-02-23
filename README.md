# isotools


command-line tools for calculating cluster quality.

## Summary

The software is intended for use with clusters formed from
of feature vectors of neural extracellular field recordings. 

These programs were tested/developed on LINUX systems, but can be compiled to
run on Microsoft Windows or Mac OS. To compile, you will need a standard C++
compiler and the make utility. Compile the C++ files from the command line with
the make command. That will produce executable files. The intermediate .o files
can be removed with make clean. make install will copy the binary files to
the top-level directory (isoitools).

## Contains
  src - source code directory, along with documentation in readme files in subdirectories
    isoi   - source code for program that calculates Isolation Information cluster quality measures
    isorat - source code for program that calculates Isolation Distance and LRatio cluster quality measures
    wave2f - source code for program that 
  wavex.txt - short recording of extracellular waveforms using a tetrode in rat CA1. The format of the
              file allows using it with wave2f. The first line of the file is the number of microwires
   	      used for recording (4, since a tetrode). The next line is the number of samples on each
	      channel (64). The following line specifies the sampling frequency (32 KHZ). The remaining
	      lines are the waveform data for each spike (4 channels at 64 samples per channel, so 256
	      numbers per line).

## Installation

To build the software, change to the src directory and use the make command. Then run the
make install command. This will copy the binaries for the 3 programs above into the top-level
directory. The usage instructions and further information on each of the programs is located
in the subdirectories as readme.txt files.

## Use

wavex.txt is a sample data file that can be used with these programs. 

The following instructions demonstrate how to use the programs:
 1. use wave2f to convert the waveforms into feature vectors & save results to file:
  wave2f wavex.txt wavefeat.txt
 2. use isoi to obtain the Isolation Information quality measures & save results to file:
  isoi wavefeat.txt wave_IsoI.txt
 3. use isorat to obtain the Isolation Distance and LRatio quality measures & save results to file:
  isorat wavefeat.txt wave_IsoDLRat.txt

Please consult the individual readme.txt files for more information on the
individual programs (located in the src  subdirectories).

## Contact

For more information contact Sam Neymotin ( samn at neurosim dot downstate dot edu )
 or Andre Fenton ( afenton at nyu dot edu ).

## References:
 The methods used in these program are described in an article in press
  at The Journal of Neuroscience,
  Measuring the quality of neuronal identification in ensemble recordings
  by Neymotin SA, Lytton WW, Olypher AO, Fenton AA (2011).




## license

<license info>
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
</license info>
