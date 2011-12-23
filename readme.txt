Author: Vlad I. Morariu
E-mail: morariu(at)umd(.)edu
Date:     2007-06-26
Modified: 2008-01-29

This is a new version of the IFGT library that is based on the original
IFGT library and FIGTree code by Vikas C. Raykar and Changjiang Yang.
Improvements include bug fixes, a new parameter selection algorithm, 
the addition of the tree-based data structure, and a C interface.

The library can be compiled both as a set MEX DLL files to provide a
Matlab interface and also as a DLL that provides a C interface.  The
original code from the IFGT and FIGTree libraries has been repackaged
into figtree.cpp with some modifications (bug fixes, extensions, and
style changes).  See figtree.h and figtree.cpp for a more detailed list of
changes.

Sample code for using the library in Matlab or C/C++ is provided in
'samples' directory.

Documentation for the Matlab interface is available in the *.m stub
files in the 'matlab' directory.  The Windows MEX files have a *.dll
extension since most Matlab versions accept *.dll files (whereas
others accept *.mex and not *.mexw32 or vice versa).  In case of
conflict with the figtree.dll that provides the C interface, the figtree.dll
MEX file can be renamed to have either a *.mex or *.mexw32 extension
(depending on the version of Matlab).

For description of available C functions, see 'include/figtree.h'.

This code makes use of the Approximate Nearest Neighbors (ANN) library
by default to provide the tree-based datastructure, but if the
appropriate flags are supplied (see below), the library can compile
without ANN.  A slightly modified version of ANN (only makefiles are
modified) is provided along with the figtree code for convenience.  For newer
releases of ANN or for any other information regarding ANN, please go to
the authors' website:

http://www.cs.umd.edu/~mount/ANN/


-----------------------------------------------------------------------
Using Precompiled WIN32 binaries
-----------------------------------------------------------------------

The figtree precompiled binaries (.dll) and lib files are provided in 
figtree/bin and figtree/lib.

The ANN precompiled binaries and lib files are provided in 
ann_1.1.1/MS_Win32/bin and ann_1.1.1/MS_Win32/lib.

To use the figtree library, you must add {BASE_PATH}/figtree/lib to the 
lib path, and add figtree.lib as a dependency.

To execute code linked against figtree (and ANN), you must make sure that 
the locations of ANN.dll and figtree.dll are in the system path.

NOTE: See below for an example of how the PATH environment variable can be
used to add the paths of libraries without having to copy them in the 
directory of each executable that uses them.


-----------------------------------------------------------------------
Using Precompiled WIN32 matlab MEX files
-----------------------------------------------------------------------

Because the MEX files are dynamically linked against ANN.dll, ANN.dll
needs to be present in the system executable path.  If ANN.dll is not
in the path, then calling any of the MEX files in MATLAB will give an
'Invalid MEX file' message.  

NOTE: using addpath() seems to work only for telling MATLAB where to find 
the MEX file, but not where to find its dependencies (in our case ANN.dll).
One way to fix this issue is to add the location of ANN.dll to the PATH
environment variable.  For example, assume the ANN.dll is located in 
C:\ann_1.1.1\MS_Win32\bin.  To add it to the PATH environment variable (this
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
      '%PATH%; C:\ann_1.1.1\MS_Win32\bin' 'in the Variable value:' field.
      The '%PATH%' is there so that the ANN.dll directory is APPENDED to the
      System variables instance of 'PATH' (which you cannot change if you are
      not an administrator).
   b) if PATH is among the user variables, select it, and then append your
      path to the list, using ';' as delimiter.  
6. If MATLAB was currently running, it must be restarted for the PATH
   variable change to take effect (until then you will still get an 
   'Invalid MEX file' error.


-----------------------------------------------------------------------
Compiling Library and Matlab MEX files in Windows using VS8
-----------------------------------------------------------------------

You can compile the library using the Visual Studio Solution
(FIGTREE_DIR)/vs8/figtree.sln.  Compiling the entire solution will compile both the
C library (figtree.dll will be in the bin directory, and figtree.lib will be in
the lib directory), and the mex files, which will be placed in the 
matlab directory.

To compile the library, the ANN library will be needed.
You must set the paths so that the *.h, *.lib, and *.dll are found 
during the compilation process.

-----------------------------------------------------------------------
Compiling Library in Linux or Solaris
-----------------------------------------------------------------------

1. (optional) Compile ANN library using g++ (see ANN compile instructions).  
   This Will generally require the command 'make linux-g++' or 
   'make sunos5-g++' for a static library and 'make linux-g++-sl' or
   'make sunos5-g++-sl' for shared libraries.  If there are errors about
   unexpected end of line in the makefile and you are using Solaris, 
   try 'gmake' instead.  

   The option 'linux-g++-sl' does not exist in version 1.1.1 of ANN but
   can be added with minimal work.  The version of ANN provided with the
   figtree code includes this change. See section 'Adding linux shared library
   target to ANN Makefile' below for details on how to modify the original
   version of the ANN source code if you wish to do so.

2. Compile FIGTree library as a shared library using the command 'make' 
   if you performed step 1 and compiled ANN as a shared library, or
   'make FIGTREE_NO_ANN=true' if you did not perform step 1 or wish to 
   compile the FIGTree library without ANN support (the 'tree' part of 
   FIGTree will not function).

   If you instead want to compile the FIGTree library as a static library, 
   then the command is 'make FIGTREE_LIB_TYPE=static' if you also compiled 
   ANN as a static library in step 1.  If you do not wish to use the ANN 
   library, then the command is 
   'make FIGTREE_LIB_TYPE=static FIGTREE_NO_ANN=true'.

Note: if you compile ANN as a static library but FIGTree as a shared library
   you might have 'relocation' errors from the linker since the static library
   is not compiled as position independent code (-fpic or -fPIC flags).  If 
   linking succeeds, then the resulting code might be slow.

3. If shared libraries are created, make sure that they are
   in the library path.  Assuming that both ANN and FIGTree libraries are 
   in your home directory, one way to do this is to set LD_LIBRARY_PATH:
    
   setenv LD_LIBRARY_PATH ${HOME}/ann_1.1.1/lib:${HOME}/figtree/lib

   if LD_LIBRARY_PATH was not previously set, or

   setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH:${HOME}/ann_1.1.1/lib:${HOME}/figtree/lib

   if LD_LIBRARY_PATH was previously set.

-----------------------------------------------------------------------
Adding linux shared library target to ANN Makefile
-----------------------------------------------------------------------

1. in ann_1.1.1/Make-config, add the following lines in the config
   option list (be careful with tabs as they have special functions
   in makefiles -- you can simply copy/paste following lines):

#					Linux using g++, make shared library
linux-g++-sl:
	$(MAKE) targets \
	"ANNLIB = libANN.so" \
	"C++ = g++" \
	"CFLAGS = -fpic -O3" \
	"MAKELIB = g++ -shared -o" \
	"RANLIB = true"

2. in ann_1.1.1/Makefile, after the comment 'main make entry point'
   add, "linux-g++-sl" in the list of possible targets.  For example
   if the list is 'alpha-g++ macosx-g++ linux-g++' then make it
   'alpha-g++ macosx-g++ linux-g++ linux-g++-sl'.

3. (optional) in ann_1.1.1/Makefile, add the following line under the 'default:'
   target after the line containing '@echo "  make linux-g++ ...':

	@echo "  make linux-g++-sl         for Linux and g++, make shared libs"


-----------------------------------------------------------------------
Compiling Matlab MEX files in Linux or Solaris
-----------------------------------------------------------------------

This has been tested on Matlab 7 with g++ as a compiler.

1. If ANN support is desired, compile ANN library as a shared library as
   described above.  It is important that the library is compiled as 
   a shared library, since code might otherwise not link or run slowly.

2. Set up matlab to use the gcc compiler by running 'mex -setup'.

3. in Matlab, set the working directory to (FIGTREE_DIR)/matlab/ and run
   CompileMexFilesUsingMatlab(use_ann, ann_header_dir, ann_lib_dir).
   If no arguments are passed, use_ann is assumed to be 1 (so step 1 must 
   be performed, and the ANN header and library directories are assumed to be
   the default locations. To compile the MEX files without ANN support you 
   can use command 'CompileMexFilesUsingMatlab(0)'.  See 
   CompileMexFilesUsingMatlab.m for more details.
