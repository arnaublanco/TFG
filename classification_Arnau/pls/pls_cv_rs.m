function RMSECV=pls_cv_rs(X,Y,prep,Pt,Nrep)
% 
% function RMSECV=pls_cv_rs(X,Y,prep,Pt,Nrep)
% 
% The function performs Random selection cross-validation
%
% INPUT ARGUMENTS
% ================
% X --> X-block. Matrix of responses. Matrix of features for the independent variables (m). Size n x m.
% Y --> Y-block. Matrix of concentrations. Matrix of features for the dependent variables (p). Size n x p.
% prep --> Preprocessing step to be applied to X and Y blocks within cross-validation steps:    
%          0: No preprocessing (DEFAULT)
%          1: Mean substraction (columns)
%          2: Autoscaling to unit variance (columns)
% Pt --> Percentage of samples which will be selected randomly for training
% from the original matrix X. Pt is set to 70% by default
% Nrep --> Number of repetitions to perform random selection. Nrep is set
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

if size(X,1)~=size(Y,1)
    error('The number of samples (rows) must be equal in X-block and Y-block');    
end

if nargin < 5
    Nrep=10;
    if nargin < 4
        Pt=70;        
        if nargin<3
            prep=0; %no preprocessing by default
        end
    end
end


%Initialization
    X0=X;
    Y0=Y;

    [N,M]=size(X);

%---------------------------------------------------------------------
%Random selection
%---------------------------------------------------------------------
Pt=Pt/100;
Nt=floor(N*Pt);

%Maximum number of LV to construct the PLS model
    max_pls=min([Nt,M]);
    if max_pls > 20 %Maximum number of LV to perform CV    
        max_pls=20;        
    end
    

%Initialization
    RMSECV=zeros(Nrep,max_pls);

for ival=1:Nrep

disp(['Random selection CV ',num2str(ival),'/',num2str(Nrep)]);
    
    %Random vector
        irnd=randperm(N);
        
    %Data for training
        X=X0(irnd(1:Nt),:);
        Y=Y0(irnd(1:Nt),:);
    
    %Data for validation    
        Xval=X0(irnd(Nt+1:end),:);    
        Yval=Y0(irnd(Nt+1:end),:);    
        Yval0=Yval;

    %Size of training data
        [nx,m]=size(X);
        [ny,py]=size(Y);
        
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
        X=X-repmat(mux,Nt,1);    %mean substraction
        Y=Y-repmat(muy,Nt,1);    %mean substraction
        Xval=Xval-repmat(mux,size(Xval,1),1);    %mean substraction
%        Yval=Yval-repmat(muy,size(Yval,1),1);    %mean substraction   
        if prep>1
            X=X./repmat(sx,Nt,1); %scale to unit variance
            Y=Y./repmat(sy,Nt,1); %scale to unit variance
            Xval=Xval./repmat(sx,size(Xval,1),1); %scale to unit variance
%            Yval=Yval./repmat(sy,size(Yval,1),1); %scale to unit variance       
        end
    end


    %Basic NIPALS algorithm [1]
    %------------------------------
    %  X=T*P'+E
    %  Y=U*Q'+F

    [T0,W0,P0,U0,B0,Q0]=pls_nipals(X,Y,max_pls);

    %=======================================================================
    %Aqui hay que sobreescribir las matrices U, Q, P, T, W, B segun el nº de LV
    %escogido
    % max_pls a partir de aqui sera el nº de LV escogido o sugerido por el programa
    %=======================================================================
    for Npls=1:max_pls

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
            RMSECV(ival,Npls)=sqrt(sum(aux_res(:).^2)/numel(aux_res));

    end

end
%---------------------------------------------------------------------
%Random selection
%---------------------------------------------------------------------

