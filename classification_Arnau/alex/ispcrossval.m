function [ ValRate, sValRate, ModelRate, sModelRate, yloo] = ispcrossval(x, y, Voptions, Coptions,grouping)
% Usage  [ ValRate, sValRate, ModelRate, sModelRate]  = ispcrossval(x, y, k, funcval,kpb, wnorm)
%		Perform crossvalidation using
%			k sets the validation:
%			k=0  leave one out validation
%			k=integer k-fold cross validation
%			k='rNUM'  random subsampling
%				where NUM is the percent to take
%                       kpb is passed to funcval as model parameter
%                       funcval: string containing a _function compatible with
%                [yValPred, yModelPred] = funcval(xTrain, yTrain, ...
%						  xVal,...)
%        type_clas: Type of clasificator
%                   'knn': k nearest neighbour[ ValRate, sValRate,
%                   ModelRate, sModelRate,yloo] = ispcrossval(x, y, k, Coptions)
%                   'svm': support vector Machine
%
% Alex Perera, eSense Systems-UB, Aug2 001
%2010: Add auto-scale for train data and substract main of train data to
%val data and divide stdtrain,
%2011-Update
%By Ana Guaman
%[ ValRate, sValRate, ModelRate, sModelRate,yloo] = ispcrossval(x, y, k,
%Coptions)
%x:data set [M samples x N features]
%y:labels [M samples x C number of clases] ej. 2 samples and 2 clases-->
%[0 1; 1 0];
%k:Type of crossvalidation: same as before
%Coptions: Structure with information about preprocessing and classifiers.
%Coptions.scale: 'Autoscale' or 'Mean Center';
%Coptions.name: Name of classifier--> 'knnval', 'fuzzyknnval', 'classify'
%Coptions.k: the number o K neighbours only knnval or fuzzyknnval
%Coptions.type: in the case of LDA the type lineal
%OUTPUT The same as before

global verbose
if isempty(verbose)
    verbose=0;
end

if nargin < 5
    grouping=0;
end

if nargin < 4
    Coptions.name='knnval';
    Coptions.k=3;
end

yloo=[];
% vecs: Number of samples
% nd: Number of dimensions
[vecs , dummy, nd] = size(x);

% nclas: Number of label classes
[yvecs, nclas] = size(y);
plott=0;

if (vecs ~= yvecs)
    error('ispcrossval: The number of rows in the data (%d rows) and the labels (%d rows) are inconsistent',vecs,yvecs);
end


Voptions=ispConvertVoptions(Voptions);

%% Validate preprocessing:
if not(isfield(Coptions,'scale'))
    Coptions.scale='nothing';
end

if not(ischar(Coptions.scale))
    error('Error in Coptions.scale');
end

%
switch lower(Coptions.scale)
    case lower('autoscale')
        Poptions=struct('name','autoscale');
    case lower('meancenter')
        Poptions=struct('name','meancenter');
    case ''
        Poptions=struct('name','');
    case lower('nothing')
        Poptions=struct('name','');
    otherwise
        error('Preprocessing %s in Coptions.scale not implemented',Coptions.scale);
end

%%
ix = 1;
lastlength=0;
chars={'|','/','-','\\'};

%numv=floor(vecs/k);

ValRate = zeros(size(y,2));
ModelRate = zeros(size(y,2));
pcValData=zeros(size(y,2));
%pcModel =zeros(size(y,2));

[kit,trainindexes,validindexes]=ispGetValidationIndex(x,grouping,Voptions,y);

iic=0;
for ss=1:kit,
    iic=iic+1;
    valInd=validindexes{ss};
    trainInd=trainindexes{ss};
    
    if nd>1
        xVal = x(valInd, :,:);
    else
        xVal = x(valInd, :);
    end
    yVal = y(valInd, :);
    if nd>1
        xTrain = x(trainInd, :,:);
    else
        xTrain = x(trainInd, :);
    end
    yTrain = y(trainInd, :);
    
    
    %Poptions.mean=mean(xTrain,1);
    %Poptions.std=std(xTrain,0,1);
    
    [xtrainauto, Poptions] = isppreprocessing(xTrain,Poptions);
    xvalauto   = isppreprocessing(xVal, Poptions);
    
    if plott
        figure();
        subplot(1,2,1),
        plot(valInd, ones(size(valInd)),'r+');hold on
        plot(trainInd,ones(size(trainInd)),'bo');
        subplot(1,2,2),
        plot(xVal(:,1), xVal(:,2),'r+');hold on
        plot(xTrain(:,1),xTrain(:,2),'bo');
        pause
    end;
    
    [yValPred,pcModel,dummy]=Dclassifiers(xtrainauto,yTrain,xvalauto,Coptions);
    
    %Calculate the Rates
    %      if (yStd2yRow(yVal) ~= yValPred) fprintf('\n'); randp(ss),fprintf('\n');
    %      end
    if strcmpi(Voptions.name,'LeaveOneOut') && nargout > 4
        yloo(validindexes{ss})=yValPred;
    end
    
    pcValData = pcValData + succ (yVal, yValPred);
    
    
    if verbose > 0
        ix = mod(ix,3)+1;
        strtp=[chars{ix} ' (' num2str(ss) '/' num2str(kit) '):' num2str(classrate(pcValData,1))];
        fprintf(repmat('\b',1,lastlength));
        fprintf(strtp);
        lastlength=length(strtp);
        %fprintf('ispcrossval: %d/%d\t:',ss,kit);
    end
    
end

if verbose > 0
    fprintf('\n');
end
ValRate = pcValData;
sValRate = pcValData;
ModelRate = pcModel;
sModelRate = pcModel;
