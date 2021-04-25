close all;
clear all;

load betas.mat
load gp.mat
load gp_test.mat

[coeff, score, ~, ~, ~, mu] = pca(train,'Centered',true);

figure;
title('Train');
scatter3(score(:,1),score(:,2),score(:,3),[],gp,'filled');

figure;
title('Train');
scatter3(score(:,1),score(:,2),score(:,3),[],gp,'filled');

% scoreTest95 = (test2-mu)*coeff(:,1:3);
% 
% figure;
% scatter3(score(:,1),score(:,2),score(:,3),[],gp,'filled');
% hold on
% scatter3(scoreTest95(:,1),scoreTest95(:,2),scoreTest95(:,3),[],gp_test,'*');