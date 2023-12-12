function [maxVar, maxInd, vars] = maxVarMap(Map)
% input is a 3D matrix (2D maps concatenated in the third dimention)

for ind=1: size(Map,3)
    map=Map(:,:,ind);
    vars(ind)=var(map(:));
end

[maxVar, maxInd]=max(vars);