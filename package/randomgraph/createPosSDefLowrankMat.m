function [T,U] = createPosSDefLowrankMat(n,r,zero_threshold)
% createPosSDefLowrankMat - creates a simple d-regular undirected graph

eigvals=[1,sort((rand(1,r-2)*(1-zero_threshold^2)+zero_threshold^2).^(1/2),'descend'),zero_threshold];

[U,~,~]=svds(randn(n,n),r);

T = (U.*eigvals)*U';

end


