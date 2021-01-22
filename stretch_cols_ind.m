% Function that normalizes a vector column by column.
%  INPUT:
%   · input: Vector to normalize.
%   · m1: Lower boundary.
%   · m2: Upper boundary.
%  OUTPUT:
%   · out: Normalized vector.
%   · pars: Minimum and maximum values of ´input´.

function [out, pars] = stretch_cols_ind(input, m1, m2)

[nRows, nCols] = size(input);
new = zeros(nRows,nCols); 
pars = zeros(nCols,2);  %% max and min for each rescaling

for i = 1:nCols
    [new(:,i), pars(i,:)] = stretch(input(:,i),m1,m2);
end

out = new;