% Function that rescales an array with a given minimum and maximum.
%  INPUT:
%   · input: Input array.
%   · m1: Minimum (lower bound)
%   · m2: Maximium (upper bound)
%  OUTPUT:
%   · out: Input array rescaled.

function [out,pars] = stretch_cols_ind(input,m1,m2)

[nRows, nCols] = size(input);
new = zeros(nRows,nCols); 
pars = zeros(nCols,2);
for i = 1:nCols
    [new(:,i), pars(i,:)] = stretch(input(:,i),m1,m2);
end

out = new;