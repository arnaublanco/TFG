function RMSECV=pls_cv(X,Y,mode,prep,Pt,Nrep)
% 
% function RMSECV=pls_cv(X,Y,mode,prep,[Pt or Kfold],[Nrep or nothing])
% 
% INPUT ARGUMENTS
% ================
% X --> X-block. Matrix of responses. Matrix of features for the independent variables (m). Size n x m.
% Y --> Y-block. Matrix of concentrations. Matrix of features for the dependent variables (p). Size n x p.
% mode --> Cross-validation mode:
%      1) Leave-one-out (DEFAULT)
%      2) Random selection
%      3) K-fold
% prep --> Preprocessing step to be applied to X and Y blocks within cross-validation steps:    
%          0: No preprocessing (DEFAULT)
%          1: Mean substraction (columns)
%          2: Autoscaling to unit variance (columns)
% 
% Pt --> ((2)Random selection). Percentage of samples which will be selected randomly for training
% from the original matrix X. Pt is set to 70% by default.
%
% Kfold --> ((3)Kfold cross-validation). Number of splits of the dataset to
% perform a K-fold cross-validation. Kfold = 2 by default.
%
% Nrep --> ((2)Random selection). Number of repetitions to perform random selection. Nrep is set
% to 10 by default.
% 
% OUTPUT ARGUMENTS
% ================
% RMSECV --> Root-mean-squared error of cross validation. Size n x npls,
% where npls is the maximum possible number of latent variables to
% construct the PLS model. npls corresponds to the minimum between n and m,
% i.e min([n,m]); where n,m are the dimensions of X-block.
% Each row of RMSE corresponds to one sample and each column corresponds to
% one LV.
%
import classification.pls.*;

if nargin < 2
    error('Not enough input arguments');
end

if nargin < 4
    prep=0;%No preprocessing
    if nargin < 3
        mode=1;     %Leave-one-out by default        
    end
end

%Random selection
%-----------------
if mode==2
    if nargin< 6
        Nrep=10;
        if nargin < 5
            Pt=70;
        end
    end
end

%K-fold
%-----------
if mode==3
    if nargin<5
        Kfold=2;
    else
        Kfold=Pt;
    end
end


%Original data
X0=X;
Y0=Y;

if size(X0,1)~=size(Y0,1)
    error('The number of samples (rows) must be equal in X-block and Y-block');    
end

%Perform cross-validation
%--------------------------
switch mode
    case 1,
        RMSECV=pls_cv_loo(X,Y,prep);
    case 2,
        RMSECV=pls_cv_rs(X,Y,prep,Pt,Nrep);
    case 3,
        RMSECV=pls_cv_kfold(X,Y,prep,Kfold);
    otherwise
        warning('Cross-validation mode must be between 1-3. No cross-validation performed');
end

