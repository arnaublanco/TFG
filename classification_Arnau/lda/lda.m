function [loadings scores v ncomp]=lda(xdata,xlabel,ncomp)
%
%function [loadings scores v ncomp]=lda(xdata,xlabel,ncomp)
%
%INPUT ARGUMENTS
%=================
% xdata: data set (nsample x nfeatures)
% xlabel: labels of data set (nsample x nclass) [0 1; 1 0] (ex)
%         The number of columns is equal to the number of classes
% ncomp: Number of components to perform LDA
%
%OUTPUT ARGUMENTS
%=================
% loadings: loadings of lda
% scores : scores of lda
% v: eigenvalues from SVD decomposition
% ncomp: Number of components used to construct LDA model
%
% See also pcalda_train  pcalda_pred
%

%Data Dimensions
    [samples features]=size(xdata);
    [samples, class]=size(xlabel);

    if nargin < 3
        ncomp=class-1;
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

if ncomp~=0
    
    %Mean of all data set
    if samples>1
        mdataset=mean(xdata);
    else
        mdataset=xdata;
    end

    %Mean of each class
    meaclass=zeros(class,features);
    sbdiff=zeros(class,features);
    swdif=zeros(samples,features);

    for i=1:class
        classi=find(xlabel(:,i)==1); % find the elements of class i
        bdatai=xdata(classi,:); % data set of class i
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
    % scores=xdata*loadings;
    loadings=s;
    scores=xdata*s(:,1:ncomp);
    v=diag(v).^2;

else
%     error('ncomp=0, maybe we do not have enough samples!!');
    loadings=0;
    scores=0;
    v=0;
    ncomp=0;
end
    