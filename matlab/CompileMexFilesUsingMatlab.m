function CompileMexFilesUsingMatlab( use_ann, ann_header_dir, ann_lib_dir )
% Compiles the Mex files for using the FIGTree in Matlab.  If no inputs
% are supplied, default values are assumed (see below).
%
% Before running this file, run 'mex -setup' and make sure that gcc or g++
% are used for compiling mex files (code was only tested using gcc/g++).
%
% If the ANN library is compiled as a static library (libANN.a), then it likely
% did not have -fpic or -fPIC as flags for the compiler, which can cause
% errors when linking the static library into a mex file (which is a 
% dynamic library).  To solve this problem, there are a couple of options:
%
%    1) compile the ANN library into a shared library (libANN.so) and 
%       make sure that it is in the library path or that LD_LIBRARY_PATH
%       is set accordingly. (PREFERRED OPTION)
%
%    2) modify (ANNBASEDIRECTORY)/Make-config to ensure that -fpic or -fPIC
%       flags are used even when compiling into a static library
%
%    3) use 'ld -G -z textoff' to link to libANN.a 
%       (Matlab uses 'g++ -shared') by default when using g++ compiler, which
%       automatically sets '-z text' ld flag which does not allow for relocation)
%       Can result in slow code. 
%
% Inputs:
%   use_ann - if nonzero, then the mex files will be compiled with ANN
%             support, and will be linked against libANN.a or libANN.so.
%             By default, mex files are compiled with ANN support.
%   ann_header_dir - the directory in which ANN headers are located.
%             By default ann_header_dir = '../../ann_1.1.1/include'
%   ann_lib_dir - the directory in which the ANN library is located.
%             By default, ann_lib_dir = '../../ann_1.1.1/lib'
%
% Created:  09-12-07 by Vlad Morariu
% Modified: 10-05-07

if( ~exist('use_ann') || use_ann ~= 0 )
    use_ann = 1;
else
    use_ann = 0;
end;

if( ~isunix )
    fprintf('To compile mex files in windows, you can use the VS8 projects provided.\n');
end;

if( use_ann == 0 )
    fprintf('Compiling mex files without ANN support\n');
    mex('-v','-O','-I../include','-DFIGTREE_NO_ANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtree.cpp','-output','figtree');
    mex('-v','-O','-I../include','-DFIGTREE_NO_ANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeChooseParametersNonUniform.cpp','-output','figtreeChooseParametersNonUniform');
    mex('-v','-O','-I../include','-DFIGTREE_NO_ANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeChooseParametersUniform.cpp','-output','figtreeChooseParametersUniform');
    mex('-v','-O','-I../include','-DFIGTREE_NO_ANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeChooseTruncationNumber.cpp','-output','figtreeChooseTruncationNumber');
    mex('-v','-O','-I../include','-DFIGTREE_NO_ANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeEvaluateDirect.cpp','-output','figtreeEvaluateDirect');
    mex('-v','-O','-I../include','-DFIGTREE_NO_ANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeEvaluateDirectTree.cpp','-output','figtreeEvaluateDirectTree');
    mex('-v','-O','-I../include','-DFIGTREE_NO_ANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeEvaluateIfgt.cpp','-output','figtreeEvaluateIfgt');
    mex('-v','-O','-I../include','-DFIGTREE_NO_ANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeEvaluateIfgtTree.cpp','-output','figtreeEvaluateIfgtTree');
    mex('-v','-O','-I../include','-DFIGTREE_NO_ANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeKCenterClustering.cpp','-output','figtreeKCenterClustering');
else
    fprintf('Compiling mex files with ANN support\n');
    if( ~exist('ann_header_dir'))
        ann_header_dir = '../../ann_1.1.1/include';
        fprintf('ann_header_dir not passed in as argument, assuming default value\n');
        fprintf('ann_header_dir=%s\n', ann_header_dir);
    end;
    if( ~exist('ann_lib_dir'))
        ann_lib_dir = '../../ann_1.1.1/lib';
        fprintf('ann_lib_dir not passed in as argument, assuming default value\n');
        fprintf('ann_lib_dir=%s\n', ann_lib_dir);
    end;

    mex('-v','-O','-I../include',sprintf('-I%s',ann_header_dir),sprintf('-L%s',ann_lib_dir),'-lANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtree.cpp','-output','figtree');
    mex('-v','-O','-I../include',sprintf('-I%s',ann_header_dir),sprintf('-L%s',ann_lib_dir),'-lANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeChooseParametersNonUniform.cpp','-output','figtreeChooseParametersNonUniform');
    mex('-v','-O','-I../include',sprintf('-I%s',ann_header_dir),sprintf('-L%s',ann_lib_dir),'-lANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeChooseParametersUniform.cpp','-output','figtreeChooseParametersUniform');
    mex('-v','-O','-I../include',sprintf('-I%s',ann_header_dir),sprintf('-L%s',ann_lib_dir),'-lANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeChooseTruncationNumber.cpp','-output','figtreeChooseTruncationNumber');
    mex('-v','-O','-I../include',sprintf('-I%s',ann_header_dir),sprintf('-L%s',ann_lib_dir),'-lANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeEvaluateDirect.cpp','-output','figtreeEvaluateDirect');
    mex('-v','-O','-I../include',sprintf('-I%s',ann_header_dir),sprintf('-L%s',ann_lib_dir),'-lANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeEvaluateDirectTree.cpp','-output','figtreeEvaluateDirectTree');
    mex('-v','-O','-I../include',sprintf('-I%s',ann_header_dir),sprintf('-L%s',ann_lib_dir),'-lANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeEvaluateIfgt.cpp','-output','figtreeEvaluateIfgt');
    mex('-v','-O','-I../include',sprintf('-I%s',ann_header_dir),sprintf('-L%s',ann_lib_dir),'-lANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeEvaluateIfgtTree.cpp','-output','figtreeEvaluateIfgtTree');
    mex('-v','-O','-I../include',sprintf('-I%s',ann_header_dir),sprintf('-L%s',ann_lib_dir),'-lANN','-DFIGTREE_USE_MATLAB_MEX','../src/figtree.cpp','../src/KCenterClustering.cpp','../src/mex/mexFigtreeKCenterClustering.cpp','-output','figtreeKCenterClustering');
end;