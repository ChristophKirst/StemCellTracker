Tracker
=============

Version Jan 22 2014 
by Christoph Kirst, 
Rockefeller University, 
ckirst@rockefeller.edu


Usage
-----

The main code is runTracker.m. 

It reads all spread sheet data files from an input directory, 
tracks cells over time, and writes augmented spread sheet data files 
to an output directory. Each cells is assigned a 'Trajectory ID' that 
is consistent between times. The output files also record for each line 
(ie cell) the previous and next cell ID's that cell is mapped to.

The code assumes that the data file names have the string (case insensitive) 
'frame=number' that defines the time sequence, and a file extension of .cvs.  
See loadEmbryoDataFile.m for the column types assumed for the input data.

The usual MATLAB help command will bring up instructions for all our functions.

For more control over various internal parameters see setParameter.m


Example
-------

    param = setParmeter();
    runTracker(dirin, dirout, param)

see the sections of /Test/testEmbryoTracker.m for more examples and tests.



Algorithm
---------
    
traking is done by matching objects in subsequent frames using
a cost matrix formalism and then determining the trajectories.


