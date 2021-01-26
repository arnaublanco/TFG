function mout=pcalda_train(model)
% 
% function mout=pcalda_train(model)
% 
% Perform PCA+LDA using SVD. train = T�P'+ E
%
% Training data matrix is decomposed in scores (T) and loadings (P), using this
% scores and loadings we construct the PCA model with "npc" components. 
%
% INPUT ARGUMENTS
% ===============
% model --> Input data Structure containing the following fields:
%           .Xt --> MxN training matrix (M trials/samples, N dimensions). Each row
%                     contains a different spectrum
%           .npc --> Number of Principal Components. 2 by default.
%           .it --> Labels to differentiate among classes/substances. E.g: 3 samples -
%               3 classes:
%                   [1 0 0;... %class 1
%                    0 1 0;... %class 2
%                    0 0 1];   %class 3
%               Each row corresponds to a different sample in the training dataset, each column to a
%               different class. If it is not specified, LDA will not be
%               performed and LDA outputs arguments will be left at zero values.
%           .prep --> Preprocessing options:
%                   0: No preprocessing (DEFAULT)
%                   1: Mean substraction (columns)
%                   2: Autoscaling to unit variance (columns)
%           .cl --> Confidence level for outlier detection. (95% by default)
%         
% OUTPUT ARGUMENTS
% ================
% mout --> Output data structure containing the following fields:
%              .T --> PCA Training Scores. Each row is a different
%              score/sample in the new basis (M x npc)
%              .P --> PCA Training Loadings. Each column is a different PC
%              (N x npc)
%              .Tlda --> LDA Training Scores. Each row is a different score/sample in the new basis
%              .Plda --> LDA Training Loadings. Each column is a different PC
%              .var --> Mx1 PCA variance vector (a column vector) in %
%              .prep --> Preprocessing applied to training data
%              .mu --> a row vector containing the mean value per each
%                     column of the training dataset "train" (1 x N)
%              .s --> a row vector containing the standard deviation per each
%                    column of the training dataset "train" (1 x N)
%              .it --> Labels of training (see input argument)
%              .pow --> PCA Recovered power by T*P' in %
%              .RMSEC --> PCA Root mean squared error between train dataset and
%                       T*P'
%              .powlda --> LDA Recovered power by Tlda*Plda' in %
%              .RMSElda --> LDA Root mean squared error between train dataset and
%                       Tlda*Plda'
%
% See also pcalda_pred
% 

if nargin < 1
    error('Not enough input arguments');
end

%Load data from input model
    train=model.Xt;
if ~isfield(model,'npc')        
    npc=2;
else    
    npc=model.npc;
end
    it=model.it;
if ~isfield(model,'prep')    
    prep=1;
else
    prep=model.prep;
end
%Confidence level for outlier detection
if ~isfield(model,'cl')
    cl=95;
else
    cl=model.cl;
end

[M,N] = size(train);
mu=mean(train,1);
s=std(train,0,1);
Nsus=size(it,2);

if npc > min([M,N])
    npc=min([M,N]);
end

% %=====================================================================
% % Suggested Number of components by CROSS-VALIDATION
% %=====================================================================
% 
% %Leave-one-out cross-validation
% %-------------------------------
% npc_loo=0;
% RMSECV_loo=pca_cv(model,1);
% npc_loo=find(RMSECV_loo==min(RMSECV_loo),1);
% if isempty(npc_loo)
%     npc_loo=0;
% end
% 
% %Random selection
% %------------------
% npc_rs=0;
% RMSECV_rs=pca_cv(model,2);
% npc_rs=find(RMSECV_rs==min(RMSECV_rs),1);
% if isempty(npc_rs)
%     npc_rs=0;
% end
% 
% %K-fold cross validation
% %------------------------
% npc_kfold=0;
% RMSECV_kfold=pca_cv(model,3,10);
% npc_kfold=find(RMSECV_kfold==min(RMSECV_kfold),1);
% if isempty(npc_kfold)
%     npc_kfold=0;
% end
% 
% %=====================================================================
% % CROSS-VALIDATION
% %=====================================================================


%Preprocessing
%--------------
if prep>0    
    train=train-repmat(mu,M,1);    %mean substraction
    if prep>1
        train=train./repmat(s,M,1); %scale to unit variance
    end
end
        
%PCA model
%------------------------------------------------------------
% construct the matrix Y
    if M>1
        Y = train / sqrt(M-1);
    else
        Y=train;
    end
    
% SVD does it all
    Y(isnan(Y) | isinf(Y))=1e-9;    %Avoid NaN and Inf because it doesn't work under such conditions
    [U,S,PC] = svd(Y,'econ'); %economy size
    %[U,S,PC] = svd(Y); %test
    
%Make first loading vector positive    
    PC(:,1) = PC(:,1)*sign(sum(PC(:,1)));
% calculate the variances
    S = diag(S);    
    v1=S.^2;

%     if M>1
%         v1 = S.^2/(M-1); %training eigenvalues
%     else
%         v1 = S.^2;
%     end
%-----------------------------------------------------

%Which is the optimum number of PC?
%===================================
    %N� PC suggested by the program
    v2=v1/sum(v1)*100;
    thr=1; %threshold to select NPC (%)
    aux=abs(diff(v2));
    npc_sug=find(aux(:,1)<thr,1);
%     npc_sug=min([find(min(v2)==v2),npc_sug]);
    
    %Any eigenvalue below 1 is not significant.
    %An eigenvalue of 1 means that it provides the equivalent of 1 variable
    %to the model
    
    if isempty(npc_sug)
        npc_sug=find(aux==min(aux));
    end
    

%=========================================================================
%IMPORTANT COMMENT
%Suggested number of components should be based on cross-validation graphs.
%=========================================================================

%NPC selected by the user or suggested by the program?

T = train*PC(:,1:npc); %cada fila es un score diferente
P = PC(:,1:npc); %new basis in PCA space
E = train-T*P';   %matrix of residuals

P0=sum(sum(train.^2)); %original power
RMSEC=sqrt(sum(sum(E.^2))/M/N); 
PE=sum(sum(E.^2)); %power of the matrix of residuals
pow=(P0-PE)/P0*100;   %recovered power

%Residuals Q
%------------
Q=sum(E.^2,2);
% Qcon=train*(eye(size(train,2))-PC(:,1:npc)*PC(:,1:npc)');
% Qcon2=PC(:,npc+1:size(PC,2)).^2*v1(npc+1:size(PC,2))*(M-1);

%Confidence limits for Q
% Q statistics is much more significant than T2 Hotelling, specially when data is clustered
%--------------------------------------------------------------------------

if npc >= length(v1)
    Qalfa=NaN;
else
    alfa = (2*cl-100)/100;
    O1=sum(v1(npc+1:size(PC,2)));
    O2=sum(v1(npc+1:size(PC,2)).^2);
    O3=sum(v1(npc+1:size(PC,2)).^3);
    if O1 == 0
        Qalfa=NaN;
    else
        h0=1-2*O1*O3/(3*O2^2);
        if h0<0.001
            h0=0.001;
        end
        Calfa = sqrt(2)*erfinv(alfa);
        %Outlier detection
        Qalfa=O1*(Calfa*sqrt(2*O2*h0^2)/O1 + 1 + O2*h0*(h0-1)/O1^2) ^(1/h0);
    end
end

% Q-statistics for a new preprocessed sample will be calculated as:
% Qval = xval*(eye(size(xval,2))-P*P')*xval'

%Hotelling's T^2 statistics
%---------------------------
T2=T*diag(1./v1(1:npc))*T';
% T2con=T*diag(1./sqrt(v1(1:npc)))*P';    

HT2=diag(T2);

if npc >= M
    Halfa = NaN;
else
    %Outlier detection
    alpha = (100-cl)/100;
    Halfa = npc*(M-1)/(M-npc)*ftest(alpha,npc,M-npc);
end

% T^2 Hotellings-statistics for a new preprocessed sample will be calculated as:
% HT2val = xval*P*diag(1./v1(1:npc))*P'*xval'

%------------------------------------
%Perform LDA
%------------------------------------
%By default, the number of components of LDA is (Nclass-1)
    ncomp=Nsus-1;
    if ncomp>npc
        ncomp=npc; %It doesn't make sense that Nlda is bigger than NPC
    end

if ~isempty(it) && size(it,1)==M
    [P0lda Tlda vlda Nlda]=lda(T,it,ncomp);
    
    if Nlda~=0
        Plda=P0lda(:,1:Nlda);
        Elda=(T-Tlda*Plda');

        P0b=sum(sum(T.^2)); %original power
        RMSElda=sqrt(sum(sum(Elda.^2))/M/N); 
        PElda=sum(sum(Elda.^2)); %power of the matrix of residuals
        powlda=(P0b-PElda)/P0b*100;   %recovered power
    else
        Nlda=0;
        Tlda=0;
        Plda=0;
        powlda=0;
        RMSElda=0;
    end
        

%     %Residuals Q
%     %-----------
%     Qlda=sum(Elda.^2,2);
%     
%     %Confidence limits for Q
%     % Q statistics is much more significant than T2 Hotelling, specially when data is clustered
%     %--------------------------------------------------------------------------
% 
%     if Nlda >= length(vlda)
%         Qalfa_lda=NaN;
%     else
%         alfa = (2*cl-100)/100;
%         O1=sum(vlda(Nlda+1:size(P0lda,2)));
%         O2=sum(vlda(Nlda+1:size(P0lda,2)).^2);
%         O3=sum(vlda(Nlda+1:size(P0lda,2)).^3);
%         if O1 == 0
%             Qalfa_lda=NaN;
%         else
%             h0=1-2*O1*O3/(3*O2^2);
%             if h0<0.001
%                 h0=0.001;
%             end
%             Calfa = sqrt(2)*erfinv(alfa);
%             %Outlier detection
%             Qalfa_lda=O1*(Calfa*sqrt(2*O2*h0^2)/O1 + 1 + O2*h0*(h0-1)/O1^2) ^(1/h0);
%         end
%     end
% 
%     % Q-statistics for a new preprocessed sample will be calculated as:
%     % Qval = xval*(eye(size(xval,2))-P*P')*xval'
% 
%     %Hotelling's T^2 statistics
%     %---------------------------
%     T2_lda=Tlda*diag(1./vlda(1:Nlda))*Tlda';
%     % T2con=T*diag(1./sqrt(v1(1:npc)))*P';    
% 
%     HT2_lda=diag(T2_lda);
% 
%     if Nlda >= M
%         Halfa_lda = NaN;
%     else
%         %Outlier detection
%         alpha = (100-cl)/100;
%         Halfa_lda = Nlda*(M-1)/(M-Nlda)*ftest(alpha,Nlda,M-Nlda);
%     end    
    
else
%    warning('LDA was not perform because training classes were not correctly specified');
    Nlda=0;
    Tlda=0;
    Plda=0;
    powlda=0;
    RMSElda=0;
end
%------------------------------------
%Perform LDA
%------------------------------------

%Obtain parameters for classifiers
%---------------------------------
if Nlda~=0
    centroide=zeros(Nsus,Nlda);
    sigma=zeros(Nsus,Nlda);
    W=repmat(eye(Nlda),[1 1 Nsus]);
    pW=repmat(eye(Nlda),[1 1 Nsus]);
    pW_out=repmat(eye(Nlda),[Nsus 1]);    
else %Only PCA
    centroide=zeros(Nsus,npc);
    sigma=zeros(Nsus,npc);
    W=repmat(eye(npc),[1 1 Nsus]);
    pW=repmat(eye(npc),[1 1 Nsus]);    
    pW_out=repmat(eye(npc),[Nsus 1]);        
end
    
distM(1:Nsus)=0;
for k=1:Nsus
    iit= (it(:,k)==1);
    
    if Nlda~=0
        aux=Tlda(iit,:);
    else %Only PCA has been calculated
        aux=T(iit,:);
    end
        
    centroide(k,:)=mean(aux,1);    
    sigma(k,:)=std(aux,1);
    
    %dist(k)=sqrt(sum(sigma(k,:).^2));
    if size(aux,1)>1
        %W(:,:,k)=pinv(cov(aux-repmat(centroide(k,:),size(aux,1),1)));
        W(:,:,k)=cov(aux);
        pW(:,:,k)=pinv(cov(aux));    %Es equivalente a centrar cada coordenada
    end
    distM(k)=sqrt(sigma(k,:)*pW(:,:,k)*sigma(k,:)');       
    
    st=size(pW,1);
    ini=(k-1)*st+1;
    fin=k*st;
    pW_out(ini:fin,:)=pW(:,:,k);    
end

%Construct the output argument:
%--------------------------------
%PCA
    mout.npc=npc;
    mout.T=T;
    mout.P=P;
    mout.var=v2;
    mout.eig=v1;            %eigenvalues
    mout.prep=prep;
    mout.mu=mu;
    mout.s=s;
    mout.pow=pow;
    mout.RMSEC=RMSEC;
%     mout.RMSECV_loo=RMSECV_loo;
%     mout.RMSECV_rs=RMSECV_rs;
%     mout.RMSECV_kfold=RMSECV_kfold;    
    mout.it=it;
%LDA
    mout.Nlda=Nlda;
    mout.Tlda=Tlda;
    mout.Plda=Plda;
    mout.powlda=powlda;
    mout.RMSElda=RMSElda;
%Number of suggested PC
    mout.npc_sug=npc_sug;
%     mout.npc_loo=npc_loo;
%     mout.npc_rs=npc_rs;
%     mout.npc_kfold=npc_kfold;
%Centroides for each class
    mout.centroide=centroide;
%Covariance matrix for each class
%     mout.cov=W;                           %esto no lo usa en prediccion
    mout.pcov=pW_out;
%Dispersion within classes
%     mout.sigma=sigma;                     %esto no lo usa en predicion
%Radius of the classes
    %mout.dist=dist;
    mout.distM=distM;
%Confidence level in outlier identification
    mout.cl=cl;
%Residuals
    mout.Q=Q; 
    mout.Qalfa=Qalfa;
%    mout.Qlda=Qlda;       
%     mout.Qalfa_lda=Qalfa_lda;
%Hotelling's T^2 statistics
    mout.HT2=HT2;
    mout.Halfa=Halfa;
%     mout.HT2_lda=HT2_lda;    
%     mout.Halfa_lda=Halfa_lda;    
%Outliers detection using Q residuals and T2 statistics
%     mout.Qout=find(Q>Qalfa);
%     mout.Tout=find(HT2>Halfa);
%     mout.QTout=find(Q>Qalfa & HT2>Halfa);
    mout.Qout=double(Q>Qalfa);
    mout.Tout=double(HT2>Halfa);
    mout.outlier=double(Q>Qalfa & HT2>Halfa);
    
       
    