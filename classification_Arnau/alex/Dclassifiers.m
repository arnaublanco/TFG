% Classifiers to be called for ispcrossval
function [YvalPred,pcModel,outputoptions]=Dclassifiers(xtrain,ytrain,xval,options)

if nargin == 2 && ischar(ytrain) && strcmpi(ytrain,'supported')
    YvalPred=true;
    switch lower(xtrain)
        case 'classify'
        case 'knnval'
        case 'fuzzyknnval'
        case 'svmda'
        otherwise
            YvalPred=false;
    end
    return
end

if nargin < 3
    error('The number of input parameters have to be more than 3');
end

if nargin < 4
    options.name='classify';
end
outputoptions=struct;

switch(lower(options.name))
    case lower('classify')
        ytrain1=zeros(size(ytrain,1),1);
        for i=1:size(ytrain,2)
            ytrain1(ytrain(:,i)==1)=i;
        end
        [class error1 post] =classify(xval,xtrain,ytrain1,options.type);
        %     YvalPred=zeros(size(xval,1),size(ytrain,2));
        %     for i=1:size(ytrain,2)
        %         YvalPred((class==i),i)=1;
        %     end
        YvalPred=class';
        outputoptions.type=options.name;
        outputoptions.error=error1;
        outputoptions.posteriori=post;
        pcModel=0;
    case lower('knnval')
        if not(isfield(options,'Moptions'))
            options.Moptions=ispmetrics('options');
        end
        pcModel=zeros(size(ytrain,2));
        [YvalPred YModelPred] = knnval(xtrain,ytrain,xval,options.k,options.Moptions);
        outputoptions.type=options.name;
        pcModel = pcModel + succ ( ytrain, YModelPred);
    case lower('fuzzyknnval')
        pcModel=zeros(size(ytrain,2));
        [YvalPred YModelPred] =fuzzyknnval(xtrain,ytrain,xval,options.k);
        outputoptions.type=options.name;
        pcModel = pcModel + succ ( ytrain, YModelPred);
    case lower('svmda')
        YtrainS=zeros(size(ytrain,1),1);
        parfor i=1:size(ytrain,1)
            if ytrain(i,1)==1
                YtrainS(i,1)=1;
            else
                YtrainS(i,1)=2;
            end
        end
        modelSVM=svmda(xtrain,YtrainS,options.OSVM);
        predSVM=svmda(xval,modelSVM,options.OSVM);
        YvalPred=predSVM.pred{1,2}';
        % [YvalPred,alpha,b,yValPredan]=SVM2class(xtrain,YtrainS,xval);
        % pcModel = 0;
        %                for i=1:size(YvalPred2,2)
        %                    yVal=YvalPred2(i);
        %                 if sign(yVal)==1
        %                         YvalPred(i,1:2)=1;
        %                     else
        %                         YvalPred(i,1:2)=2;
        %                 end
        %                end
        %YvalPred(sign(YvalPred)==-1)=2;
        outputoptions.type=options.name;
        outputoptions.modelSVM=modelSVM;
        outputoptions.predSVM=predSVM;
        pcModel=0;
    otherwise
        error('Unrecognized classifier name: %s"',options.name);
end

end