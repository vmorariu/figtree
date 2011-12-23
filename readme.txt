Author: Vlad I. Morariu
E-mail: morariu(at)umd(.)edu
Date:     2007-06-26
Modified: 2008-12-05

FIGTree is a library that provides a C/C++ and MATLAB interface
for speeding up the computation of the Gauss Transform.  Its  
authors are Vlad Morariu, Vikas Raykar, and Changjiang Yang.  The 
authors worked under the supervision of Professor Ramani Duraiswami 
and Professor Larry Davis, at the University of Maryland.

The current version of the library is based on the following paper:

Vlad I. Morariu, Balaji Vasan Srinivasan, Vikas C. Raykar, 
Ramani Duraiswami, and Larry S. Davis. Automatic online tuning for 
fast Gaussian summation. Advances in Neural Information Processing 
Systems (NIPS), 2008.

For more up-to-date information on this library, please go to the 
FIGTree homepage: 

http://www.umiacs.umd.edu/~morariu/figtree/ 

or the SourceForge download page: 

http://sourceforge.net/projects/figtree

This code makes use of the Approximate Nearest Neighbors (ANN) library
for the tree-based datastructure, but if the appropriate flags are
supplied (see below), the library can compile without ANN.  A slightly
modified version of ANN (only makefiles are modified) is provided along 
with the figtree code for convenience.  For newer releases of ANN or for 
any other information regarding ANN, please go to the authors' website:

http://www.cs.umd.edu/~mount/ANN/


-----------------------------------------------------------------------

Samples and Documentation

-----------------------------------------------------------------------

Sample code for using the library in Matlab or C/C++ is provided in
'samples' directory.

Documentation for the Matlab interface is available in the *.m stub
files in the 'matlab' directory.  The Windows MEX files have a *.dll
extension since most Matlab versions accept *.dll files (whereas
others accept *.mex and not *.mexw32 or vice versa).  In case of
conflict with the figtree.dll that provides the C interface, the 
figtree.dll MEX file can be renamed to have either a *.mex or *.mexw32
extension (depending on the version of Matlab).

For a description of available C/C++ functions, see 'include/figtree.h'.


-----------------------------------------------------------------------

Using Precompiled WIN32 binaries

-----------------------------------------------------------------------

The figtree precompiled binaries (.dll) and lib files are provided in 
figtree/bin and figtree/lib.

To use the figtree library, you must add {BASE_PATH}/figtree/lib to the 
lib path, and add figtree.lib as a dependency.

To execute code (not MATLAB) linked against figtree (and ANN), you must 
make sure that the locations of ann_figtree_version.dll and figtree.dll 
are in the system path.

NOTE: See below for an example of how the PATH environment variable can 
be used to add the paths of libraries without having to copy them in
the directory of each executable that uses them.



-----------------------------------------------------------------------

Compiling Library and Matlab MEX files in Windows using VS8

-----------------------------------------------------------------------

You can compile the library using the Visual Studio Solution
(FIGTREE_DIR)/vs8/figtree.sln.  Compiling the entire solution will 
compile both the C/C++ library (figtree.dll will be in the bin 
directory, and figtree.lib will be in the lib directory), and the mex
files, which will be placed in the matlab directory.

Make sure you first switch to the 'Release' mode, or else figtree will
be very slow!



-----------------------------------------------------------------------

Compiling Matlab MEX files in Linux, Solaris, or Windows using Matlab mex

-----------------------------------------------------------------------

This has been tested on Matlab 7 with g++ as a compiler.  Also, it works
in windows if the compiler is set up right.

1. Set up matlab to use the gcc compiler by running 'mex -setup'.

2. In Matlab, set the working directory to (FIGTREE_DIR)/matlab/ and run
   CompileMexFilesUsingMatlab.m.



-----------------------------------------------------------------------

Compiling Library in Linux or Solaris

-----------------------------------------------------------------------

1. Compile FIGTree library as a shared library by issuing command 'make'
   
   If you instead want to compile the FIGTree library as a static library, 
   then the command is 'make FIGTREE_LIB_TYPE=static'.

2. If shared libraries are created, make sure that they are
   in the library path.  Assuming that FIGTree libraries are 
   in your home directory, one way to do this is to set LD_LIBRARY_PATH:
    
   setenv LD_LIBRARY_PATH ${HOME}/figtree/lib

   if LD_LIBRARY_PATH was not previously set, or

   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HOME}/figtree/lib

   if LD_LIBRARY_PATH was previously set.



-----------------------------------------------------------------------

Setting the PATH variable for programs that link against figtree.lib

-----------------------------------------------------------------------

Assume the figtree.dll is located in 
C:\figtree\bin.  To add it to the PATH environment variable (this
does not require admin privileges):

1. Right click on 'My Computer'
2. Click on 'Properties'
3. Click on the 'Advanced' tab
4. At the bottom, click on 'Environment Variables'
5. In the upper section 'User variables for (yourusername)' is a 
   list of variables with their values.  Variable 'PATH' may or may
   not be one of them.
   a) if PATH is not among the user variables (usually TEMP and TMP are 
      present), then you must add it by clicking 'New', entering 'PATH'
      in the 'Variable name:' field, and then entering 
      '%PATH%; C:\figtree\bin' 'in the Variable value:' field.
      The '%PATH%' is there so that the figtree directory is APPENDED to the
      System variables instance of 'PATH' (which you cannot change if you are
      not an administrator).
   b) if PATH is among the user variables, select it, and then append your
      path to the list, using ';' as delimiter.  




-----------------------------------------------------------------------

Todo List for Future versions

-----------------------------------------------------------------------

- TODO: Clean up makefiles for Unix/Linux/Solaris version.

- TODO: estimate avg neighbors using approach of Faloutsos et 
          al (SIGMOD 2000) instead of querying kd-tree on subsampled
          sources
- TODO: currently, k-center clustering is performed using an
          N log(K) + K^2 (expected time) approach based on the Gonzales 
          K-center clustering.  Use N log(K) approach instead.
- TODO: add function to improve cost prediction by tuning to actual 
          hardware (i.e. floating point operation times vary relative 
          to memory access times on different machines; also exp() 
          takes different number of float ops depending on machine)
- TODO: add memory allocation checks to avoid crashes
- TODO: allow kd-tree, clusters, and coefficients computed on sources 
          to persist between function calls in case user wants to reuse 
          them
- TODO: make some more of the internal parameters easier to modify by 
          users without requring code to be compiled again (such 
          parameters include: amt of subsampling to do for 
          estimating cost, truncation number limit, etc)
- TODO: modify ANN so that it no longer uses global variables (which
          prevents figtree from being
          called simultaneously in multiple threads)

-----------------------------------------------------------------------

Change Log

-----------------------------------------------------------------------

-----------------------------------------------------------------------
figtree-0.9.2  2008-12-05  (Changes made by Vlad Morariu)  
-----------------------------------------------------------------------

- NEW: Added function that chooses fastest method between Direct, 
    DirectTree, Ifgt, and IfgtTree, making figtree a black-box approach
- Added point-wise and cluster-wise approaches to selecting 
    source/target truncation numbers
- Moved ann-1.1.1 library inside (figtree-dir)/external directory to 
    reduce confusion
- Change mex code to include ANN code (instead of linking to it) to
    prevent issues with library paths
- Removed many of the mex files that were really just confusing people, 
    and left only figtree, method selection, and k-center clustering
    mex files as they are most useful.
- Changed interface of the mex files to remove redundant arguments.
- Other improvements (see source files)

-----------------------------------------------------------------------
figtree-0.9.1  2008-02-26  (Changes made by Vlad Morariu)
-----------------------------------------------------------------------

- FIX:  parameter selection gave bad params since K (number of 
        clusters) was not allowed to equal number of sources

-----------------------------------------------------------------------
figtree-0.9    2008-01-30   (Changes made by Vlad Morariu)
-----------------------------------------------------------------------

- Initial adaptation from Vikas's initial code
- NEW:  integrated Vikas's FIGTree code into library
- NEW:  added parameter selection method that does not assume uniformly 
    distributed sources
- NEW:  added C/C++ interface, with one function (figtree()) that does
    most work
- Sources and targets no longer need to fit in unit hypercube

