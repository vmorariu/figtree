% Script the demonstrate the use of FIGTree.  Note, only 
% figtree/matlab/figtree.dll is needed for the figtree() function to 
% function.  All other mex files are provided in case users require
% advanced control over parameter selection and clustering processes.  See
% sample2.m for example of how to use the additional mex files provided.
clear all;
close all;
clear functions;
clc;

% the MEX binaries for the FIGTree library reside in ../matlab
% we assume that the script is run from the figtree/samples/ directory,
% otherwise the correct path for the mex files must be set.
addpath('../matlab');

disp('---------------------------------------------');
disp(sprintf('Examples to demonstrate the use of FIGTree'));
disp('---------------------------------------------');

% the number of sources
N = 5000;
disp(sprintf('Number of source points N=%d',N));

% the number of targets
M = 5000;
disp(sprintf('Number of target points M=%d',M));

% the desired error
epsilon = 1e-3;     
disp(sprintf('Target error epsilon=%e',epsilon));

% run all algorithms with all parameter selection methods with some 
% different source/target distributions, dimensionality, and bandwidth, 
% to show some cases where some methods perform better than others
for i = 1:4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate sources and targets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear X_shifted_scaled X_shifted X Y;
    fprintf('\n\n');
    if( i == 1 )
        disp('----------------------------------------------------------------------');
        disp(sprintf('Sources drawn from multiple gaussians\nTargets uniformly distributed'));
        
        % the data dimensionality
        d = 4;
        disp(sprintf('Dimensionality d=%d',d));
       
        % the bandwidth
        h = .01; %sqrt(2)*0.2*sqrt(d);
        disp(sprintf('Bandwidth h=%f',h));
        disp('----------------------------------------------------------------------');
        
        % the source points
        % d x N matrix of N source points in d dimensions.
        % Scale the data
        G = 25; 
        m = rand(d,G);
        v = 0.02*ones(1,G);
        [X] = generate_multiple_gaussians(N,G,m,v,d);
 
        for j = 1:d
            shift = min(X(j,:));
            X_shifted(j,:) = X(j,:) - shift;
            scale = 1/max(X_shifted(j,:));
            X_shifted_scaled(j,:) = X_shifted(j,:)*scale;
        end

        X = X_shifted_scaled;

        % the target points
        % d x M matrix of M source points in d dimensions.
        Y = rand(d,M);
    end;
    
    if( i == 2 )
        disp('----------------------------------------------------------------------');
        disp(sprintf('Sources and targets drawn from (same) multiple gaussians'));
        
        % the data dimensionality
        d = 4;
        disp(sprintf('Dimensionality d=%d',d));
       
        % the bandwidth
        h = .1; %sqrt(2)*0.2*sqrt(d);
        disp(sprintf('Bandwidth h=%f',h));       
        disp('----------------------------------------------------------------------');
        
        % the source points
        % d x N matrix of N source points in d dimensions.
        % Scale the data
        G = 25; 
        m = rand(d,G);
        v = 0.02*ones(1,G);
        [X] = generate_multiple_gaussians(N+M,G,m,v,d);

        for j = 1:d
            shift = min(X(j,:));
            X_shifted(j,:) = X(j,:) - shift;
            scale = 1/max(X_shifted(j,:));
            X_shifted_scaled(j,:) = X_shifted(j,:)*scale;
        end

        inds = randperm(M+N);
        X = X_shifted_scaled(:,inds(1:N));
        Y = X_shifted_scaled(:,inds(N+1:N+M));       
    end;

    if( i == 3 )
        disp('----------------------------------------------------------------------');
        disp(sprintf('Sources and targets drawn from (same) multiple gaussians'));
        
        % the data dimensionality
        d = 2;
        disp(sprintf('Dimensionality d=%i',d));
       
        % the bandwidth
        h = .1; %sqrt(2)*0.2*sqrt(d);
        disp(sprintf('Bandwidth h=%f',h));       
        disp('----------------------------------------------------------------------');
        
        % the source points
        % d x N matrix of N source points in d dimensions.
        % Scale the data
        G = 25; 
        m = rand(d,G);
        v = 0.02*ones(1,G);
        [X] = generate_multiple_gaussians(N+M,G,m,v,d);

        for j = 1:d
            shift = min(X(j,:));
            X_shifted(j,:) = X(j,:) - shift;
            scale = 1/max(X_shifted(j,:));
            X_shifted_scaled(j,:) = X_shifted(j,:)*scale;
        end

        X = X_shifted_scaled(:,1:N);
        Y = X_shifted_scaled(:,N+1:N+M);       
    end;

    if( i == 4 )
        disp('----------------------------------------------------------------------');
        disp(sprintf('Sources and targets drawn from (same) multiple gaussians'));
        
        % the data dimensionality
        d = 2;
        disp(sprintf('Dimensionality d=%i',d));
       
        % the bandwidth
        h = .2; %sqrt(2)*0.2*sqrt(d);
        disp(sprintf('Bandwidth h=%f',h));       
        disp('----------------------------------------------------------------------');
        
        % the source points
        % d x N matrix of N source points in d dimensions.
        % Scale the data
        G = 25; 
        m = rand(d,G);
        v = 0.02*ones(1,G);
        [X] = generate_multiple_gaussians(N+M,G,m,v,d);

        for j=1:d
            shift = min(X(j,:));
            X_shifted(j,:) = X(j,:)-shift;
            scale = 1/max(X_shifted(j,:));
            X_shifted_scaled(j,:) = X_shifted(j,:)*scale;
        end

        X = X_shifted_scaled(:,1:N);
        Y = X_shifted_scaled(:,N+1:N+M);       
    end; 
    
    % the source weights
    % 1 x N row vector
    q = rand(1,N);  

    % scale and shift the data and bandwidth to make sure that such 
    %   changes do not affect accuracy or speed.  (Thus, sources
    %   and targets no longer have to lie in a unit hypercube).
    new_scale = .01;
    new_shift = 0;
    X = new_scale*X+new_shift;
    Y = new_scale*Y+new_shift;
    h = new_scale*h;

    % try for various scales
    for j = 1:9


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluate gauss transform
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % direct method (exact up to machine precision)
        tic; 
        [ g_direct           ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 0, 0 ); 
        t_direct(j) = toc;

        % direct evaluation with approximate nearest neighbor structure 
        % (parameter selection not relevant, so it can be anything) - 
        % implemented as in the FIGTree paper (Raykar et al, submitted to 
        % Joural of Machine Learning).
        tic;
        [ g_direct_tree       ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 2, 0 );
        t_direct_tree(j) = toc;

        % truncated series, uniform parameter selection - implemented exactly
        % as described in IFGT paper (Raykar et al, 2005).
        tic; 
        [ g_ifgt_u      ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 1, 0 ); 
        t_ifgt_u(j) = toc;

        % truncated series, with nonuniform parameter selection - the
        % default method when calling function without last two arguments
        % that select the method: figtree( d, N, M, 1, X, h, q, Y, epsilon );
        % works well in general
        tic; 
        [ g_ifgt_nu     ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 1, 1 );
        t_ifgt_nu(j) = toc;

        % truncated series with approximate nearest neighbor structure,
        % with uniform parameter selection - implemented as in the FIGTree
        % paper (Raykar et al, submitted to Joural of Machine Learning).
        tic;
        [ g_ifgt_tree_u  ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 3, 0 );
        t_ifgt_tree_u(j) = toc;

        % truncated series, with non-uniform parameter selection.  Is faster 
        % than ifgt-nu if targets are not identically distributed to 
        % sources.
        tic;
        [ g_ifgt_tree_nu ] = figtree( d, N, M, 1, X, h, q, Y, epsilon, 3, 1 );
        t_ifgt_tree_nu(j) = toc;

        % show errors
        e_direct_tree(j)       = max(abs(g_direct-g_direct_tree))/sum(q);
        e_ifgt_u(j)      = max(abs(g_direct-g_ifgt_u))/sum(q);
        e_ifgt_nu(j)     = max(abs(g_direct-g_ifgt_nu))/sum(q);
        e_ifgt_tree_u(j)  = max(abs(g_direct-g_ifgt_tree_u))/sum(q);
        e_ifgt_tree_nu(j) = max(abs(g_direct-g_ifgt_tree_nu))/sum(q);
        
        new_scale = 4;
        new_shift = 10^j;
        X = new_scale*X+new_shift;
        Y = new_scale*Y+new_shift;
        h = new_scale*h;
    end;
    
    % show times and errors
    fprintf('direct       : %3.2e seconds\n', t_direct );
    fprintf('direct-tree  : %3.2e seconds, speedup = %3.2f\n', [t_direct_tree; t_direct./t_direct_tree] );
    fprintf('ifgt-u       : %3.2e seconds, speedup = %3.2f\n', [t_ifgt_u; t_direct./t_ifgt_u] );
    fprintf('ifgt-nu      : %3.2e seconds, speedup = %3.2f\n', [t_ifgt_nu; t_direct./t_ifgt_nu] );
    fprintf('ifgt-tree-u  : %3.2e seconds, speedup = %3.2f\n', [t_ifgt_tree_u; t_direct./t_ifgt_tree_u] );
    fprintf('ifgt-tree-nu : %3.2e seconds, speedup = %3.2f\n', [t_ifgt_tree_nu; t_direct./t_ifgt_tree_nu] );
    fprintf('\nMaximum absolute error from direct (exact) method:\n');
    fprintf('direct-tree  : %3.2e\n', e_direct_tree );
    fprintf('ifgt-u       : %3.2e\n', e_ifgt_u );
    fprintf('ifgt-nu      : %3.2e\n', e_ifgt_nu );
    fprintf('ifgt-tree-u  : %3.2e\n', e_ifgt_tree_u );
    fprintf('ifgt-tree-nu : %3.2e\n', e_ifgt_tree_nu );
end;
