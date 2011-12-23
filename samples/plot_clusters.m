function plot_clusters(N,d,X,K,ClusterIndex,ClusterCenter)
% Plots the output of k-center clustering

K=double(K);
ClusterIndex=double(ClusterIndex);

if d ==2
    sat=[0.1:0.8/K:0.9];
    %figure;
    
    for i=1:K
        ind=find(ClusterIndex==(i-1));
        col_hsv=[sat(i) sat(i) 1];
        h=plot(X(1,ind),X(2,ind),'.');
        set(h,'MarkerEdgeColor',hsv2rgb(col_hsv));
       % set(h,'Markersize',6);
       % set(h,'Linewidth',3);
        hold on;
        h=plot(ClusterCenter(1,i),ClusterCenter(2,i),'kx');
        set(h,'Markersize',8);
        set(h,'Linewidth',2);       
    end
    
    %Plot the cluster centers
    
    axis([-0.1 1.1 -0.1 1.1]);
    box on;
elseif d==3
    sat=[0.1:0.8/K:0.9];
    
    for i=1:K
        ind=find(ClusterIndex==(i-1));
        col_hsv=[sat(i) sat(i) 1];
        h=plot3(X(1,ind),X(2,ind),X(3,ind),'.');
        set(h,'MarkerEdgeColor',hsv2rgb(col_hsv));
        set(h,'Markersize',2);
        set(h,'Linewidth',2);
        hold on;
        h=plot3(ClusterCenter(1,i),ClusterCenter(2,i),ClusterCenter(3,i),'kx');
        set(h,'Markersize',8);
        set(h,'Linewidth',2);       
    end
    box on;
    
    
else
    disp('Plotting supported only in two and three dimensions.');
end
