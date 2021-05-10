% Function that computes ANOVA per voxel.
%  INPUT:
%   · data: fMRI data.
%   · gp: Labels.
%  OUTPUT:
%   · res: Results from ANOVA.

function [res] = voi_ANOVA(data, gp)

mSize = size(data,2);
res = zeros(mSize+1,2);   % one extra for ANOVA on mean betas

% Compute one way ANOVA per voxel
for i = 1:mSize
    [res(i,1), anovatab] = anova1(data(:,i),gp,'off');  % p-val
    res(i,2) = anovatab{2,5}; %f-val
end

% Compute the ANOVA on the average betas
mdata = mean(data,2);  %% row means (over voxels)
[res(mSize+1,1), anovatab] = anova1(mdata,gp,'off');
res(mSize+1,2) = anovatab{2,5};

    

