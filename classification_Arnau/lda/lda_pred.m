function [predictions, extra_info]=lda_pred(x, model)
%
%function [model]=lda_pred(x,labels,options)
%
% INPUT ARGUMENTS
% =================
% x: data set to predict a class from (nsample x nfeatures)
% model: LDA model trained with lda_train
%
% OUTPUT ARGUMENTS
% =================
% predictions: Predicted probabilities
% extra_info: Structure with:
%  - scores: Scores of LDA for the given samples
%
% AUTHORS
% ========
% 2014, Sara Rica
% 2014, Sergio Oller

import classification.lda.*;

if nargin == 0
    help lda_pred
    return
end

if nargin < 2
    error('Model is missing');
end

nsamplespred = size(x,1);

scores = x * model.loadings(:,1:model.ncomp);

dmah=zeros(nsamplespred, model.nclass);
for k=1:model.nclass
    dmah(:,k)=disteusq(model.centroide(k,:),scores,'xs',model.pcov{k}')';
end

%Prediction: Assign class to each sample
%Outlier detection
%------------------
%What samples are within K sigma?
%Samples outside K sigma will be classified as NONE (all 0)

dmah2 = dmah;
dmah2(dmah>model.ksigma*repmat(model.distM',[nsamplespred,1])) = Inf;

% We can use the softmax transform to convert the distances to the
% centroids into "probabilities"

expdmah = exp(-dmah2);
predictions = expdmah ./ repmat(sum(expdmah,2),[1,model.nclass]);


% [mdmah,imdmah]=min(dmah2,[],2);
% 
% labels=zeros(nsamplespred, model.nclass); %class assigned by the classifier
% for k=1:nsamplespred
%     if ~isinf(mdmah(k))
%         labels(k,imdmah(k))=1;
%     end
% end

extra_info = struct();
extra_info.scores = scores;
extra_info.dist_mah = dmah;
extra_info.dist_mah_outlier = model.ksigma*model.distM';

