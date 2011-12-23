% Script the demonstrate the use of FIGTREE without the wrapper function 'figtree'
% (i.e. by manually calling figtreeChooseParameters*(), KCenterClustering(),
% etc).  Manually calling these functions instead of letting figtree() do it
% allows the user to change parameters (such as number of clusters,
% truncation number, or even the clustering method itself).
clear all;
close all;
clear functions;
clc;

% the MEX binaries for the FIGTREE library reside in ../matlab
% we assume that the script is run from the figtree/samples/ directory,
% otherwise the correct path for the mex files must be set.
addpath('../matlab');

disp('---------------------------------------------');
disp(sprintf('Example to demonstrate the use of FIGTree'));
disp('---------------------------------------------');

% the data dimensionality
d = 2;
disp(sprintf('Dimensionality d=%d\n',d));

% the number of sources
N = 5000;
disp(sprintf('Number of source points N=%d\n',N));

% the number of targets
M = 5000;
disp(sprintf('Number of target points M=%d\n',M));

% the source points
% d x N matrix of N source points in d dimensions.
% Scale the data
G = 25; 
m = rand(d,G);
v = 0.02*ones(1,G);
[X] = generate_multiple_gaussians(N,G,m,v,d);

for j = 1:d
    shift = min(X(j,:));
    X_shifted(j,:) = X(j,:)-shift;
    scale = 1/max(X_shifted(j,:));
    X_shifted_scaled(j,:) = X_shifted(j,:)*scale;
end
X= X_shifted_scaled;

% the target points
% d x M matrix of M source points in d dimensions.
Y = rand(d,M);

% the source weights
% 1 x N row vector
q=rand(1,N);  

% the bandwidth
h = sqrt(2)*0.2*sqrt(d);
disp(sprintf('Bandwidth h=%f\n',h));

% the desired error
epsilon=1e-3;     
disp(sprintf('Target error epsilon=%e\n',epsilon));

% The upper limit on the number of clusters.
kLimit = N; % we'll let KCenterClustering stop whenever it reaches N clusters
            % or when there are no more unique source pts
disp(sprintf('kLimit=%d\n',kLimit));
disp(sprintf('Press any key to continue...\n'));
pause

disp('---------------------------------------------');
disp(sprintf('Choosing the FIGTree parameters\n'));
disp('---------------------------------------------');

% Choose the parameters
%
% K     -- number of clusters
% p_max -- maximum truncation number
% r     -- cutoff radius

to=clock;
% you can uncomment the next line and comment the one after it to choose
% parameters assuming uniform source distribution
%[K,p_max,r]=figtreeChooseParametersUniform(d,h,epsilon,kLimit, 1);
[K,p_max,r]=figtreeChooseParametersNonUniform(d,N,X,h,epsilon,kLimit,1);
parameters_time=etime(clock,to);

disp(sprintf('Number of clusters K=%d\n',K));
disp(sprintf('Maximum truncation number p_max=%d\n', p_max));
disp(sprintf('Cutoff radius r=%f\n',r));
disp(sprintf('Time taken=%f secs\n',parameters_time));
disp(sprintf('Press any key to continue...\n'));
pause

disp('---------------------------------------------');
disp(sprintf('Running the k-center clustering\n'));
disp('---------------------------------------------');

% k-center clustering
to=clock;
[K, rx, clusterIndex, clusterCenter, numPoints, clusterRadii] = figtreeKCenterClustering(d, N, X, K);
clustering_time=etime(clock,to);

disp(sprintf('Actual clusters(K)=%f\n', K));
disp(sprintf('Maximum cluster radius=%f\n', rx));
disp(sprintf('Time taken=%f secs\n', clustering_time));

%Pretty plot of the results of the clustering procedure. 
%Plots only for two and three dimensions.
plot_clusters(N, d, X, K, clusterIndex, clusterCenter);

disp(sprintf('Press any key to continue...\n'));
pause

disp('---------------------------------------------');
disp(sprintf('Updating the truncation number\n'));
disp('---------------------------------------------');

to = clock;
[p_max] = figtreeChooseTruncationNumber(d, h, epsilon, rx, 1);
trunc_time = etime(clock,to);

disp(sprintf('Updated Maximum Truncation Number=%d\n',p_max));
disp(sprintf('Time taken=%f secs\n',trunc_time));
disp(sprintf('Press any key to continue...\n'));
pause

disp('---------------------------------------------');
disp(sprintf('Running the  FIGTree\n'));
disp('---------------------------------------------');

to=clock;
% you can use another evaluation method by uncommenting it and commenting
% the current one.  Note that to use the last, the bandwidth should be set
% to some very small number (eg. .01) or else evaluation will take a very
% long time.
[G_FIGTREE]=figtreeEvaluateIfgt(d, N, M, 1, X, h, q, Y, p_max, K, ...
                                clusterIndex, clusterCenter, clusterRadii, r, epsilon);
%[G_FIGTREE]=figtreeEvaluateIfgtTree(d, N, M, 1, X, h, q, Y, p_max, K, ...
%                                    clusterIndex, clusterCenter, clusterRadii, r, epsilon);
%[G_FIGTREE]=figtreeEvaluateDirectTree(d, N, M, X, h, q, Y, epsilon);  % very slow if bandwidth is above h = .1
FIGTREE_time=etime(clock,to);

disp(sprintf('Time taken=%f secs\n',FIGTREE_time));
disp(sprintf('Press any key to continue...\n'));

pause

disp('---------------------------------------------');
disp(sprintf('Running the direct method.\n'));
disp('---------------------------------------------');

to=clock;
[G_direct]=figtreeEvaluateDirect(d,N,M,X,h,q,Y);
GT_time=etime(clock,to);

disp(sprintf('Time taken=%f secs\n',GT_time));

disp('---------------------------------------------');
disp(sprintf('Summary\n'));
disp('---------------------------------------------');

FIGTREE_total_time=parameters_time+clustering_time+trunc_time+FIGTREE_time;
FIGTREE_err=max(abs((G_direct-G_FIGTREE)))/sum(q);

disp(sprintf('Direct computation takes %f secs\n',GT_time));
disp('---------------------------------------------');
disp(sprintf('FIGTree takes %f secs Speedup=%f\n',FIGTREE_total_time,GT_time/FIGTREE_total_time));
disp(sprintf('Actual error for FIGTree is %e Target was %e \n',FIGTREE_err,epsilon));
disp('---------------------------------------------');



