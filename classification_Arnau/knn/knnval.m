function [yValPred] = knnval(xTrain, yTrain, xVal, k, Moptions)
% Usage [yValPred] = knnval(xTrain, yTrain, xVal, k,Moptions)
%     This Function trains a Set using xTrain yTrain
%         and estimates xVal using the resulting model
%       Using knn Classifier
%        if k is not given, a k=3 default index is given
% INPUT:
% ======
%   xTrain: Matrix nsamples_train x nfeatures.
%   yTrain: Matrix of binary labels (nsamples_train x nclasses)
%   xVal: Matrix nsamples_prediction x nfeatures
%   k: Number of neighbours (default: 3)
%   Moptions: Metric options
%
% OUTPUT:
% =======
%   yValPred: Predicted binary labels nsamples_prediction x nclasses
%
% See also: ispmetrics, ispclassify_train, ispclassify_pred, knn1
%
% AUTHORS:
% 2003, Alex Perera
% 2014, Sara Rica
% 2014, Sergio Oller
import classification.knn.*;
if nargin < 4
    k=3;
end

if nargin < 5
    Moptions=ispmetrics('options');
end

% To be estimated somehow..

[yValPred] = knn1(xTrain, yTrain, k, xVal, Moptions);
end