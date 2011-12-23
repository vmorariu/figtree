% Script the demonstrate the use of FIGTree.  Note, only 
% (FIGTREE-DIRECTORY)/matlab/figtree.dll is needed for the figtree() function to 
% function. 
clear all;
close all;
clear functions;
clc;

% the MEX binaries for the FIGTree library reside in ../matlab
% we assume that the script is run from the figtree/samples/ directory,
% otherwise the correct path for the mex files must be set.
addpath('../matlab');

% This example runs all evaluation methods with both parameter selection 
% methods with some different source/target distributions, dimensionality,
% and bandwidth, to show some cases where some methods perform better
% than others.
%
% As expected, the non-uniform parameter selection outperforms the 
% uniform parameter selection method as it tunes parameters to the actual
% data.  In addition, depending on the situation, different evaluation
% methods outperform, but the automatic method selection approach is always
% very close to the fastest method. 
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate gauss transform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % direct method (exact up to machine precision) (evalMethod=0)
    tic; 
    [ g_direct           ] = figtree( X, h, q, Y, epsilon, 0 ); 
    t_direct = toc;
    fprintf('direct       : %3.2e seconds\n', t_direct );
   
    % direct evaluation with tree data structure (evalMethod=2)
    tic;
    [ g_direct_tree       ] = figtree( X, h, q, Y, epsilon, 2 );
    t_direct_tree = toc;
    fprintf('direct-tree  : %3.2e seconds, speedup = %3.2f\n', t_direct_tree, t_direct/t_direct_tree );
      
    % ifgt (evalMethod=1) and parameter selection that assume
    % uniform distribution (ifgtParamMethod=0)
    tic; 
    [ g_ifgt_u      ] = figtree( X, h, q, Y, epsilon, 1, 0 ); 
    t_ifgt_u = toc;
    fprintf('ifgt-u       : %3.2e seconds, speedup = %3.2f\n', t_ifgt_u, t_direct/t_ifgt_u );

    % ifgt (evalMethod=1) and parameter selection that does not assume
    % uniform distribution (ifgtParamMethod=1)
    tic; 
    [ g_ifgt_nu     ] = figtree( X, h, q, Y, epsilon, 1, 1 );
    t_ifgt_nu = toc;
    fprintf('ifgt-nu      : %3.2e seconds, speedup = %3.2f\n', t_ifgt_nu, t_direct/t_ifgt_nu );

    % ifgt with tree data structure (evalMethod=3) and parameter selection 
    % that assumes uniform distribution (ifgtParamMethod=0)
    tic;
    [ g_ifgt_tree_u  ] = figtree( X, h, q, Y, epsilon, 3, 0 );
    t_ifgt_tree_u = toc;
    fprintf('ifgt-tree-u  : %3.2e seconds, speedup = %3.2f\n', t_ifgt_tree_u, t_direct/t_ifgt_tree_u );

    % ifgt with tree data structure (evalMethod=3) and parameter selection 
    % that does NOT assume uniform distribution (ifgtParamMethod=1)
    tic;
    [ g_ifgt_tree_nu ] = figtree( X, h, q, Y, epsilon, 3, 1 );
    t_ifgt_tree_nu = toc;
    fprintf('ifgt-tree-nu : %3.2e seconds, speedup = %3.2f\n', t_ifgt_tree_nu, t_direct/t_ifgt_tree_nu );

    % automatically select evaluation method (evalMethod=4) and use
    % parameter selection that assumes uniform distribution 
    % (ifgtParamMethod=0).
    tic;
    [ g_auto_u ] = figtree( X, h, q, Y, epsilon, 4, 0 );
    t_auto_u = toc;
    fprintf('auto-u       : %3.2e seconds, speedup = %3.2f\n', t_auto_u, t_direct/t_auto_u ); 
    
    % automatically select evaluation method (evalMethod=4) and use
    % parameter selection that does NOT assume uniform distribution 
    % (ifgtParamMethod=1). 
    tic;
    [ g_auto_nu ] = figtree( X, h, q, Y, epsilon, 4, 1 );
    t_auto_nu = toc;
    fprintf('auto-nu      : %3.2e seconds, speedup = %3.2f\n', t_auto_nu, t_direct/t_auto_nu );

    
    err_direct_tree = max(abs(g_direct-g_direct_tree))/sum(q);
    err_ifgt_u = max(abs(g_direct-g_ifgt_u))/sum(q);
    err_ifgt_nu = max(abs(g_direct-g_ifgt_nu))/sum(q);
    err_ifgt_tree_u = max(abs(g_direct-g_ifgt_tree_u))/sum(q);
    err_ifgt_tree_nu = max(abs(g_direct-g_ifgt_tree_nu))/sum(q);
    err_auto_u = max(abs(g_direct-g_auto_u))/sum(q);
    err_auto_nu = max(abs(g_direct-g_auto_nu))/sum(q);
        
    % show errors
    fprintf('\nMaximum absolute error from direct (exact) method:\n');
    fprintf('direct-tree  : %3.2e\n', err_direct_tree);
    fprintf('ifgt-u       : %3.2e\n', err_ifgt_u);
    fprintf('ifgt-nu      : %3.2e\n', err_ifgt_nu);
    fprintf('ifgt-tree-u  : %3.2e\n', err_ifgt_tree_u);
    fprintf('ifgt-tree-nu : %3.2e\n', err_ifgt_tree_nu);
    fprintf('auto-u       : %3.2e\n', err_auto_u);
    fprintf('auto-nu      : %3.2e\n', err_auto_nu);
    
    if( (err_direct_tree < epsilon) && ...
        (err_ifgt_u < epsilon) && ...
        (err_ifgt_nu < epsilon) && ...
        (err_ifgt_tree_u < epsilon) && ...
        (err_ifgt_tree_nu < epsilon)&& ...
        (err_auto_u < epsilon)&& ...
        (err_auto_nu < epsilon) )
        fprintf('\nDesired error (%e) satisfied.\n', epsilon);
    else
        fprintf('\nDesired error (%e) NOT satisfied!\n', epsilon);
    end;
    
end;