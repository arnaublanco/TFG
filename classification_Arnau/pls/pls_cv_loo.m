function RMSE=pls_cv_loo(X,Y,prep)
% 
% function RMSE=pls_cv_loo(X,Y,prep)
% 
% The function performs leave-one-out cross-validation
%
% INPUT ARGUMENTS
% ================
% X --> X-block. Matrix of responses. Matrix of features for the independent variables (m). Size n x m.
% Y --> Y-block. Matrix of concentrations. Matrix of features for the dependent variables (p). Size n x p.
% prep --> Preprocessing step to be applied to X and Y blocks within cross-validation steps:    
%          0: No preprocessing (DEFAULT)
%          1: Mean substraction (columns)
%          2: Autoscaling to unit variance (columns)
% 
% OUTPUT ARGUMENTS
% ================
% RMSE --> Root-mean-squared error of cross validation. Size n x npls,
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

if nargin<3
    prep=0; %no preprocessing by default
end

if size(X,1)~=size(Y,1)
    error('The number of samples (rows) must be equal in X-block and Y-block');    
end

%Initialization
X0=X;
Y0=Y;
ncomp=min([size(X,1)-1,size(X,2)]);
if ncomp>20
    ncomp=20;
end
RMSE=zeros(size(X,1),ncomp);

%---------------------------------
%Leave-one-out cross-validation
%---------------------------------
for ival=1:size(X0,1)
    
disp(['LOO CV ',num2str(ival),'/',num2str(size(X0,1))]);
    
    X=X0; Y=Y0; %Original data without preprocessing
    
    %Select validation samples without preprocessing
        Xval=X(ival,:);
    %    Xval0=Xval;    
        Yval=Y(ival,:);
        Yval0=Yval;    

    %Select training samples without preprocessing   
        X(ival,:)=[];
        Y(ival,:)=[];

    %Size of training data
        [nx,m]=size(X);
        [ny,py]=size(Y);
    
    %Maximum number of LV to construct the PLS model
        ncomp=min([nx,m]);
        if ncomp > 20 %Maximum number of LV to perform CV
            ncomp=20;
        end

    %Mean and standard deviations of training data
    if prep>0
        mux=mean(X,1);
        sx=std(X,0,1);
        muy=mean(Y,1);
        sy=std(Y,0,1);
    end
    
    %Preprocessing
    %--------------
    if prep>0    
        X=X-repmat(mux,nx,1);    %mean substraction
        Y=Y-repmat(muy,ny,1);    %mean substraction
        Xval=Xval-repmat(mux,size(Xval,1),1);    %mean substraction
%        Yval=Yval-repmat(muy,size(Yval,1),1);    %mean substraction   
        if prep>1
            X=X./repmat(sx,nx,1); %scale to unit variance
            Y=Y./repmat(sy,ny,1); %scale to unit variance
            Xval=Xval./repmat(sx,size(Xval,1),1); %scale to unit variance
%            Yval=Yval./repmat(sy,size(Yval,1),1); %scale to unit variance       
        end
    end


    %Basic NIPALS algorithm [1]
    %------------------------------
    %  X=T*P'+E
    %  Y=U*Q'+F

    [T0,W0,P0,U0,B0,Q0]=pls_nipals(X,Y,ncomp);

    %=======================================================================
    %Aqui hay que sobreescribir las matrices U, Q, P, T, W, B segun el nº de LV
    %escogido
    % ncomp a partir de aqui sera el nº de LV escogido o sugerido por el programa
    %=======================================================================
    for Npls=1:ncomp

%        U=U0(:,1:Npls);
        Q=Q0(:,1:Npls);
        P=P0(:,1:Npls);
%        T=T0(:,1:Npls);
        W=W0(:,1:Npls);
        B=B0(1:Npls);

    %matrix of regression coefficients
        reg=zeros(Npls*py,m);
        C = zeros(py*Npls,m);
        if py == 1
            C = (W*pinv(P'*W)*diag(B))';
        else
            Cp = (W*pinv(P'*W)*diag(B))';
            for k = 1:Npls
                Cpp = Cp(k,:);
                C((k-1)*py+1:k*py,:) = diag(Q(:,k))*Cpp(ones(py,1),:);
            end
        end
        reg(1:Npls*py,:)=C;

        if py>1
            for k=2:Npls
                jj       = (k-1)*py+1;
                i0       = jj-py;
                reg(jj:k*py,:) = reg(jj:k*py,:) + reg(i0:(k-1)*py,:);
            end
        else
            reg  = cumsum(reg,1);
        end

        reg=reg((1:py)+(py*(Npls-1)),:)';

        %Prediction
            Ypred=Xval*reg;

        %Undo preprocessing applied to Y-block
            if prep==1
                Ypred=Ypred+repmat(muy,size(Ypred,1),1);
            elseif prep==2
                Ypred=Ypred.*repmat(sy,size(Ypred,1),1)+repmat(muy,size(Ypred,1),1);
            end

        aux_res=Ypred-Yval0;
        
        %Root-mean-squared error of cross-validation
            RMSE(ival,Npls)=sqrt(sum(aux_res(:).^2)/numel(aux_res));

    end

end
%---------------------------------
%Leave-one-out cross-validation
%---------------------------------
