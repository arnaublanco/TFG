function [model]=lda_train(x, labels, options)
%
%function [model]=lda(x,labels,options)
%
%INPUT ARGUMENTS
%=================
% x: data set (nsample x nfeatures)
% labels: labels of data set (nsample x nclass) [0 1; 1 0] (ex)
%         The number of columns is equal to the number of classes
% options : Input Data structure containing the following fields.
%        ncomp: The number of components for performing the LDA (num
%        classes -1)
%        ksigma: Mahalanobis distance threshold in sigma units for outlier
%        detection (default: 4)
%
%OUTPUT ARGUMENTS
%=================
% model : Output Data structure containing the following fields.
%   -model.loadings: loadings of lda
%   -model.scores : scores of lda
%   -model.v: eigenvalues from SVD decomposition
%   -model.ncomp: Number of components used to construct LDA model
%   -model.centroide: Matrix of size (nclass x ncomp) with the average
%     position in the score space of the samples of each class
% 
import classification.lda.*;

if nargin == 0
    help lda_train
    return
end

if nargin == 1 && strcmpi(x, 'options')
    model = struct();
    model.ncomp = -1; %
    model.ksigma = 4; % Mahalanobis distance threshold in sigma units for outlier detection 
    return
end

if nargin < 2
    error('Labels are missing');
end

if nargin < 3
    options = lda_train('options');
end

%Data Dimensions
[samples, features]=size(x);
[samples, class]=size(labels);

if options.ncomp == -1
    ncomp=class-1; % Number of classes - 1
else
    ncomp = options.ncomp;
end

if ncomp>class-1
    ncomp=class-1;
end

%Avoid the curse of dimensionality
if ncomp>floor(samples/5)
    ncomp=floor(samples/5);
end
if ncomp > min([samples features])
    ncomp = min([samples features]);
end

if ncomp == 0
    error('ncomp=0, maybe we do not have enough samples!!');
end

   
%Mean of all data set
if samples>1
    mdataset=mean(x);
else
    mdataset=x;
end

%Mean of each class
meaclass=zeros(class,features);
sbdiff=zeros(class,features);
swdif=zeros(samples,features);

for i=1:class
    classi=find(labels(:,i)==1); % find the elements of class i
    bdatai=x(classi,:); % data set of class i
    if size(bdatai,1)>1
        meaclass(i,:)=mean(bdatai); % mean of class i
    else
        meaclass(i,:)=bdatai;
    end
    %Scatter within-class(Sw)
    swdif(classi,:)=bdatai-repmat(meaclass(i,:),length(classi),1);
    %Scater between class(Sb)
    sbdiff(i,:)=(mdataset-meaclass(i,:)).*length(classi);
end

Sb=sbdiff'*sbdiff;
Sw=swdif'*swdif;
A=pinv(Sw)*Sb;
A(isnan(A) | isinf(A))=1e-9;    %Avoid NaN and Inf because it doesn't work under such conditions
[s v d]=svd(A);
% loadings=s(:,1:ncomp);
% scores=x*loadings;
loadings=s;
scores=x*s(:,1:ncomp);
v=diag(v).^2;
centroide = zeros(class,ncomp);
sigma = zeros(class, ncomp);
pcov = cell(class,1);
distM = zeros(class,1);
for i=1:class
    classi= labels(:,i)==1; % find the elements of class i
    bdatai=scores(classi,:); % data set of class i
    sigma(i,:)=std(bdatai,1);
    if size(bdatai,1)>1
        centroide(i,:)=mean(bdatai); % mean of class i
        pcov{i} = pinv(cov(bdatai));
    else
         centroide(i,:)=bdatai;
    end
    distM(i)=sqrt(sigma(i,:)*pcov{i}*sigma(i,:)');
end


model = struct();
model.name = 'lda';
model.loadings = loadings;
model.scores = scores;
model.v = v;
model.nclass = class;
model.ncomp = ncomp;
model.centroide = centroide;
model.pcov = pcov;
model.distM = distM;
model.ksigma = options.ksigma;
end
