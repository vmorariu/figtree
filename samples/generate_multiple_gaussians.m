function [X]=generate_multiple_gaussians(N,G,m,v,d)
% Generates N samples in d dimensions  from G gaussians
%
% m -- dxG mean 
%
% v  -- 1xG variance 

n=round(N/G)*ones(1,G);
n(G)=n(G)+(N-n(1)*G);

X=[];
for i=1:G
    X=cat(2,X,m(:,i)*ones(1,n(i))+v(i)*randn(d,n(i)));
end

