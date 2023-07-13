function [T,r,L] = createRandRegWeightedGraph(n, deg, mode, varargin)
%createRandRegWeightedGraph This function creates a normalized weight
%matrix T = D^(-1/2)*W*D^(-1/2), where W comes from a random regular graph
%of degree deg, the non-zero weights are chosen according to the model
%'mode', and D is a diagonal matrix containing the degrees of the graph.
L=createRandRegGraph(n,deg);
if strcmp(mode,'Uniform_far_from_zero')
    zero_threshold = varargin{1};
    r=(rand(n*deg/2,1)*(1-zero_threshold))+zero_threshold;
end
T = Reg_Diff_Matrix(r,L); 
end

