function [yValPred, yModelPred] =fuzzyknnval(xTrain, yTrain, xVal, k)
  
  % Usage [yValPred yModelPred] =fuzzyknnVal(xTrain, yTrain, xVal, k)
  %     This Functin trains a Set using xTrain yTrain
  %         and estimates xVal using the resulting model
  %         yValPred output from xVal
  %       Using fuzzy knn Classifier
  %        if k is not given, a k=3 default index is given 
  %        m parameter is 2 
  %
  % Also see: fuzzyknncrisp, fuzzyknnval
  %  
  % March 2003 Rafa Rubio
   
   uTrain=zeros(size(yTrain,1),size(yTrain,2));   

   for i=1:size(yTrain,1)
     [vv,ii]=max(yTrain(i,:));      
     uTrain(i,ii(1))=1;
   end
  
  if nargin < 4
    k=1;
  end
  
  global silent
  [vecs dims] = size(xTrain);
  [yvecs nclas] = size(yTrain);
  plott=0;
  calculate_pred=0;
  
  % To be estimated somehow..
  
  if ~silent,
    printf('KnnVal: Using k=%d \n', k );
    
    printf('knnVal: Cal. Val Perf.\n\r');
  end

  [yValPred] = fuzzyknncrisp(xTrain, yTrain, xVal, k);

  %%% inline knn
  
  
  if calculate_pred,
    if ~silent,
      printf('\n');
      printf('knnVal: Cal. Model Perf.\n\r');
      
    end
    [yModelPred] = fuzzyknncrisp(xTrain, yTrain,xTrain,k);
    else
    if ~silent,
      printf('knnVal:Don''t Cal. Model Perf.\n\r');
    end
    yModelPred = yStd2yRow(uTrain);
  end
  
  if plott,
    clf;
    subplot(1,2,1),
    plot(xVal(:,1), xVal(:,2),'r+');hold on
    plot(xTrain(:,1),xTrain(:,2),'bo');
    subplot(1,2,2),
    plotKnn(xTrain, yTrain, k,xVal);
    title('Training')
    %figure (2);clf;plotpKnn(xTrain, yRow2 yStd(yModelPred),...
			     %   k);
    %title('Model')
    %figure (3);clf;plotpknn(xVal, yRow2yStd(yValPred), ...
			     %  k );
    %title('Validation')
    pause
  end
  
  