close all;
clear all;

load betas.mat
load betasC.mat

betas = reshape(betas,[size(betas,1) size(betas,2)*size(betas,3)]);
coeff = pca(betas);

figure;
scatter3(coeff(:,1),coeff(:,2),coeff(:,3))