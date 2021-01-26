function [labels, extra_info]=pls_pred(X,plsmodel)
% 
% function [labels, extra_info]=pls_pred(X,plsmodel)
% 
% INPUT ARGUMENTS
% ===============
% X --> X-block. Matrix of n responses. Matrix of features for the
%       independent variables (m). Size n x m.
% plsmodel --> Data structure containing at least the following fields:
%              .reg--> Matrix of regression coefficients. Size m x p. (Y=X*reg)
%              .T --> X-scores matrix. Size n x ncomp
%              .W --> X inner weights matrix. Size m x ncomp
%              .P --> X-loadings matrix. Size m x ncomp
%              .U --> Y-scores matriz. Size n x ncomp
%              .B --> Y inner regression coefficients. Size ncomp x 1
%              .Q --> Y-loadings matrix. Size p x ncomp
%              .ssq --> Residuals sum of squares (%). First column: X-block, 2nd column: Y-block
%                   Each row contains the percentage of residual for X-block or
%                   Y-block after extracting k LV.
%              .prep --> Preprocessing options applied to training data (X-Y):
%                   0: No preprocessing (DEFAULT)
%                   1: Mean substraction (columns)
%                   2: Scale to unit variance (columns)
%              .mux --> a row vector containing the mean value per each
%                     column of the training dataset "Xt"
%              .sx --> a row vector containing the standard deviation per each
%                    column of the training dataset "Xt"
%              .muy --> a row vector containing the mean value per each
%                     column of the training dataset "Yt"
%              .sy --> a row vector containing the standard deviation per each
%                    column of the training dataset "Yt"
%
% plsmodel must be obtained using "pls_train" in the training stage
%
% Note 1: The same preprocessing step as Xtrain is applied to X for
% prediction.
% Note 2: The reverse preprocessing step is applied to Y for
% prediction using the preprocessing step applied to Ytrain.
% 
% OUTPUT ARGUMENTS
% ================
% labels --> Y-block. Matrix of predicted concentrations.
%                      Matrix of features for the dependent variables (p). Size n x p.
% extra_info --> Output data structure containing different fields
%
% See also pls_train  
%
import classification.pls.*;

if nargin < 2
    error('Not enough input arguments');
end

%Load input data
    reg=plsmodel.reg;
    W=plsmodel.W;
    P=plsmodel.P;
    prep=plsmodel.prep;
    mux=plsmodel.mux;
    sx=plsmodel.sx;
    muy=plsmodel.muy;
    sy=plsmodel.sy;

ncomp=size(P,2);    
[n,m]=size(X);

%Preprocessing
%--------------
if prep>0    
    X=X-repmat(mux,n,1);    %mean substraction    
    if prep>1
        X=X./repmat(sx,n,1); %scale to unit variance
    end
end

%Prediction
    Y=X*reg;

%Project new data onto the basis specified by X-loadings P
    T=X*W*pinv(P'*W);

%Undo preprocessing applied to Y-block
if prep==1
    Y=Y+repmat(muy,n,1);
elseif prep==2
    Y=Y.*repmat(sy,n,1)+repmat(muy,n,1);
end
    
%-----------------------------------
%Outlier detection using PLS model
%-----------------------------------

    %Qresiduals
        Qval=diag(X*(eye(m)-plsmodel.P*plsmodel.P')*X');
    %T^2 Hotelling/Leverage statistics
    if n>1
        f=diag(sqrt(1./(diag(T'*T)/(n-1))));
    else
        f=diag(sqrt(1./(diag(T'*T))));
    end

    if ncomp>1
        tsq1=(sum((f*T').^2))';
    else
        tsq1=((f*T').^2)';
    end

    if n>1
        leverage=tsq1/(n-1);
    else
        leverage=tsq1;
    end

    iQval=(Qval>plsmodel.Qalfa);
    ilevval=(leverage>plsmodel.Halfa);

%    outlier = (iQval & ilevval); %No me convence esta manera de calcular outliers

%-----------------------------------
%Outlier detection using PLS model
%-----------------------------------
    
%Output data structure
    labels=Y;                   %Predicted
    extra_info=struct();
    extra_info.T=T;                   %Scores in LV space
%Outlier analysis
    extra_info.Qx=Qval;
    extra_info.leverage=leverage;    
    extra_info.Qout=double(iQval);
    extra_info.levout=double(ilevval);
    extra_info.outlier=double(iQval & ilevval);       %Outlier (No me convence como se calcula)
    