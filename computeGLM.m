% Function that computes the General Linear Model.
%  INPUT:
%   · data: Input data (voxels x timepoints).
%   · dm: Design matrix.
%   · zTimeSeries: 1 -> z-score data; 0 -> no z-scoring of data
%   · zBetas: 1 -> z-score betas; 0 -> no z-scoring of betas

function [betas,t] = computeGLM(data, dm, zTimeSeries, zBetas)

if(zTimeSeries==1)
    data = zscore(data);  % z-score data (voxel-wise)
end

nVox = size(data,2);
betas = zeros(size(dm,2),nVox);
t = zeros(size(betas));

% Independently for each voxel, fit GLM
for vox = 1:nVox
    [B,~,stats] = glmfit(dm,data(:,vox),'normal','constant','off');
    betas(:,vox) = B;
    t(:,vox) = stats.t;
end

% Z-score betas within each voxel, across trials
if(zBetas)
    betas = zscore(betas);  % z-score betas (voxel-wise)
end

end