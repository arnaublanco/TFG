function mout=pcalda_pred(val,model,K,none_thr)
% 
% function mout=pcalda_pred(val,model,K,none_thr)
% 
% Perform PCA+LDA using SVD to the validation dataset "val" using the training
% model "model" which can be obtained using the function "pcalda_train"
%
% INPUT ARGUMENTS
% ===============
% val --> MxN validation/test matrix (M trials/samples, N dimensions)
% model --> Structure containing the following fields (it can be obtained from "pcalda_train"):
%              .T --> PCA Training Scores. Each row is a different score/sample in the new basis
%              .P --> PCA Training Loadings. Each column is a different PC
%              .Tlda --> LDA Training Scores. Each row is a different score/sample in the new basis
%              .Plda --> LDA Training Loadings. Each column is a different PC
%              .var --> Nx1 PCA variance vector (a column vector) in %
%              .eig --> Nx1 PCA eigenvalues (a column vector) of the
%                       covariance matrix
%              .prep --> Preprocessing applied to training data
%                       0: No preprocessing (DEFAULT)
%                       1: Mean substraction (columns)
%                       2: Scale to unit variance (columns)
%              .mu --> a row vector containing the mean value per each
%                     column of the training dataset
%              .s --> a row vector containing the standard deviation per each
%                    column of the training dataset
%              .it --> Labels of training (see ival)
%              .pow --> PCA Recovered power by T*P' in %
%              .RMSE --> PCA Root mean squared error between train dataset and
%                       T*P'
%              .powlda --> LDA Recovered power by Tlda*Plda' in %
%              .RMSElda --> LDA Root mean squared error between train dataset and
%                       Tlda*Plda'
%
% K --> K-nearest neighbours parameters to be used in KNN classification (K=3 by default)
%
% none_thr --> threshold to decide when a sample will be classified as NONE

% 
% OUTPUT ARGUMENTS
% ================
% mout --> Prediction model which contains the following fields:
%            .Tval --> PCA Validation/test Scores. Each row is a different score/sample in the new basis
%            .Tval_lda --> PCA+LDA Validation/test Scores. Each row is a
%                         different score/sample in the new basis
%            .iEpred --> Classes assigned to validation samples by Euclidean classifier
%            .iMpred --> Classes assigned to validation samples by Mahalanobis classifier
%            .iKpred --> Classes assigned to validation samples by KNN
%                        classifier (K=3)
%
% See also pcalda_train
%

if nargin < 2
    error('Not enough input arguments');
end

if nargin < 3
    K=3;
end

%threshold to decide when a sample will be classified as NONE
if nargin < 4
    none_thr=50;    
end

%Extract information from PCA-LDA model
    npc = model.npc;
    Nlda = model.Nlda;
    if Nlda~=0
        Tt = model.Tlda;
    else %Only PCA has been performed
        Tt = model.T;
    end

if size(Tt,1)~=size(model.it,1)
    error('Number of training samples is not equal to the number of training labels');
end

%Maximum number of neighbours depending on the minimum number of samples in
%all the classes
    maxK=min(sum(model.it));
    if K>maxK
        K=maxK;
    end

[M,N] = size(val);
Nsus=size(model.it,2);      %N� of substances
%Ksigma=6;

%Preprocessing
%--------------
prep=model.prep;
if prep>0    
    val=val-repmat(model.mu,M,1);    %mean substraction
    if prep>1
        val=val./repmat(model.s,M,1); %scale to unit variance
    end
end

%Project validation samples to the new basis constructed using PCA
    Tval=val*model.P; %PCA projection
%Project validation samples to the new basis constructed using PCA+LDA
    if Nlda~=0
        Tval_lda=Tval*model.Plda; %LDA projection
    else
        Tval_lda=0; %Output argument when LDA has not been performed
    end
    
% confE=zeros(Nsus,M);
% confM=zeros(Nsus,M);
% confK=zeros(1,M);
% confK_c=zeros(1,M);
confK_c2=zeros(1,M);
% confP=zeros(1,M);
% pdf=zeros(1,Nsus);

%-----------------------------------
%Outlier detection using PCA model
%-----------------------------------
    eig=model.eig;        %eigenvalues

    %Qresiduals
        Qval=diag(val*(eye(N)-model.P*model.P')*val');
    %T^2 Hotelling statistics
        HT2val=diag(Tval*diag(1./eig(1:size(model.P,2)))*Tval');

iQval=(Qval>model.Qalfa);
iHT2val=(HT2val>model.Halfa);

outlier = (iQval & iHT2val);

%-----------------------------------
%Outlier detection using PCA model
%-----------------------------------

% 
% %Euclidean classifier
% %-----------------------
%     deucl=disteusq(model.centroide,Tval_lda,'xs');
%     
%     for m=1:M
%         for k=1:Nsus
%             if model.dist(k)==0
%                 confE(k,m)=0;
%             else
%                 confE(k,m)=calc_subs_conf(deucl(k,m),0,model.dist(k));
%             end
%         end
%     end
%           
%     %What samples are within K sigma?
%     %Samples outside K sigma will be classified as NONE (all 0)
%     deucl2=deucl.*(deucl<Ksigma*repmat(model.dist',1,M));
%     deucl2(deucl2==0)=NaN;
%     
%     [mde,imde]=min(deucl2);
% 
%     %Prediction: Assign class to each sample
%         ival_Epred=zeros(M,Nsus); %class assigned by the classifier
%         for k=1:M
%             if ~isnan(mde(k))
%                 ival_Epred(k,imde(k))=1;
%             end               
%         end
%         
% 

pcov=model.pcov;

%Mahalanobis classifier
%----------------------
dmah=zeros(Nsus,M);
if Nlda~=0
    for k=1:Nsus
        ini=(k-1)*Nlda+1;
        fin=k*Nlda;
        dmah(k,:)=disteusq(model.centroide(k,:),Tval_lda,'xs',pcov(ini:fin,:));
    end
else
    for k=1:Nsus 
        ini=(k-1)*npc+1;
        fin=k*npc;        
        dmah(k,:)=disteusq(model.centroide(k,:),Tval,'xs',pcov(ini:fin,:));
    end
end    

% 
%     for m=1:M
%         for k=1:Nsus
%             if model.distM(k)==0
%                 confM(k,m)=0;
%             else
%                 confM(k,m)=calc_subs_conf(dmah(k,m),0,model.distM(k));
%             end
%         end
%     end
% 
%     %Outlier detection
%     %------------------
%     %What samples are within K sigma?
%     %Samples outside K sigma will be classified as NONE (all 0)
%     dmah2=dmah.*(dmah<Ksigma*repmat(model.distM',1,M));
%     dmah2(dmah2==0)=NaN;
% 
%     [mdmah,imdmah]=min(dmah2);
%     
%     %Prediction: Assign class to each sample
%         ival_Mpred=zeros(M,Nsus); %class assigned by the classifier
%         for k=1:M
%             if ~isnan(mdmah(k))
%                 ival_Mpred(k,imdmah(k))=1;
%             end
%         end


%KNN classifier (K=3)
%--------------------
    if Nlda~=0
        dknn=disteusq(Tt,Tval_lda,'xs');
    else
        dknn=disteusq(Tt,Tval,'xs');
    end
    
    [ydknn,idknn]=sort(dknn);
    %We are interested in the first K neighbors
    idk=idknn(1:K,:);
        
    ival_Kpred=zeros(M,size(model.it,2)); %class assigned by the classifier

    for k=1:M %Evaluate each validation sample
        %Prediction: Assign class to each sample
        %----------------------------------------            
        lab=model.it(idk(:,k),:);        
        mm=max(sum(lab,1));        
        imm=find(mm==sum(lab,1));        
        if length(imm)==1       %no hay empate 
            col=imm;            
        else %hay empate => hay que desempatar        
            %Cogemos la distancia promedio a las clases seg�n el numero            
            %de vecinos que han empatado            
            daux=zeros(1,length(imm));            
            for k2=1:length(imm)            
                ilab = lab(:,imm(k2))==1;                
                aux=dknn(idk(ilab,k),k);                
                daux(k2)=mean(aux,1);                
            end            
            [mdaux,imdaux]=min(daux);
            col=imm(imdaux);            
        end
        %class assigned
%        if ~isnan(mdmah(k))         %If the sample is not an outlier
        if outlier(k) == 0
            ival_Kpred(k,col)=1;             
            
            %Confidence in substance identification
            %--------------------------------------
            %Distance to the nearest neighbour of the assigned class
%             lab1=find(lab(:,col)==1,1);
%             dknn2=dknn(idk(lab1,k),k);      %Esto es una distancia euclidea!!!
            if model.distM(col)~=0
%                 confK(k)=calc_subs_conf(dknn2,0,model.distM(col));  %Distancia euclideo con sigma de mahalanobis
%                 confK_c(k)=calc_subs_conf(dmah(col,k),0,model.distM(col)); %distancia de mahalanobis con sigma de mahalanobis
            %This confidence takes into account Mahalanobis dispersion.
                d1=dmah(:,k);
                d2=1./d1;
                d3=d2./repmat(sum(d2,1),Nsus,1);
            %What happens if a distance is exactly 0? (on the centroid)
                d3(isnan(d3))=1;
                confK_c2(k)=d3(col)*100;
                
%             %Probabilistic confidence
%             %========================
%             for iclass=1:Nsus
%                 pdf(iclass)=exp( -0.5*(Tval_lda(k,:)-model.centroide(iclass,:))*model.pcov(:,:,iclass)*(Tval_lda(k,:)-model.centroide(iclass,:))' )...
%                              *sqrt(abs(det(model.pcov(:,:,iclass)))) / (2*pi)^(size(model.Plda,2)/2);
%             end
%             prob=pdf/sum(pdf);
%             confP(k)=prob(col)*100;
%             %========================
            
            %Sample classified as NONE
            %==========================
            if confK_c2(k)<none_thr       
                ival_Kpred(k,col)=0;
            end
            %==========================
            
%             else
%                 if dknn2==0
%                     confK(k)=100;
%                 else
%                     confK(k)=0;
%                 end
            end
            
        else      %If the sample is an outlier
%             confK(k)=0;
%             confK_c(k)=0; %This confidence takes into account Mahalanobis dispersion.
            confK_c2(k)=0; %This confidence takes into account Mahalanobis dispersion. (1/d)
%             confP(k)=0; %Probabilistic confidence based on Mahalanobis distance
        end
        
    end
        

%Construct output argument
%--------------------------
%Scores
    mout.Tval=Tval;
    mout.Tval_lda=Tval_lda;
%Classification results        
%     mout.iEpred=ival_Epred;    
%     mout.iMpred=ival_Mpred;    
    mout.iKpred=ival_Kpred;
%Confidence in substance classification
%     mout.confE=confE';
%     mout.confM=confM';
%     mout.confK=confK;   
%     mout.confK_c=confK_c;   
    mout.conf=confK_c2;
%     mout.confP=confP;
%K parameter
    mout.K=K;
%Mahalanobis distance training data
%     mout.dmah=dmah;         %esto me lo puedo cargar
%     mout.d3=d3;             %esto me lo puedo cargar
%Residuals
    mout.Qval=Qval;
    mout.HT2val=HT2val;    
%Outliers
    mout.outlier=double(outlier);
    mout.none=none_thr;
    

