% Function that computes one-way ANOVA (univariate ANOVA) for each voxel in VOI.
% INPUT:
%   - data: Matrix of betas.
%   - gp: Array containing stimuli numbering (1: Forest, 2: People, 3:
%   Traffic).
% OUTPUT:
%   - res: Vector containing p-value and F-test statistic of ANOVA per
%   voxel.

function [res] = voi_ANOVA(data, gp)

mSize = size(data,2);
res = zeros(mSize+1,2); % One extra row for ANOVA on mean betas

% Compute one-way ANOVA per voxel
for i = 1:mSize
    [res(i,1), anovatab] = anova1(data(:,i),gp,'off'); % Compute ANOVA (´anova1´: Returns p-value and ANOVA table)
    res(i,2) = anovatab{2,5}; % F-test statistic
end

% Compute the ANOVA on the average betas
mdata = mean(data,2);  % Row mean (mean over voxels)
[res(mSize+1,1), anovatab] = anova1(mdata,gp,'off'); % Compute ANOVA
res(mSize+1,2) = anovatab{2,5}; % F-test statistic

    

