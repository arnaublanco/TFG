generatedata=0;
%options.metric='Euclidean';
options.metric='Mahalanobis';

%% Test Euclidean
if generatedata==1
    nsamplestrain=1000;
    nsamplestest=200;
    nfeats=50;
    nlabels=10;
    options.k=3;
    
    nsamples=nsamplestrain+nsamplestest;
    data=zeros(nsamples,nfeats);
    labels=zeros(nsamples,nlabels);
    labels1d=randi(nlabels,nsamples,1);
    for i=1:nlabels
        labels(labels1d==i,i)=1;
        featsfromlabel=randi(nfeats,randi(nfeats,1,1),1);
        myrandom=(rand-0.5)*2;
        data(labels1d==i,featsfromlabel)=data(labels1d==i,featsfromlabel)+myrandom*10;
        
    end
    labelspred_official=knn1(data(1:nsamplestrain,:),labels(1:nsamplestrain,:),k,data(nsamplestrain+1:nsamples,:));
    save(sprintf('testknn1_data_%s.mat',options.metric));
else
    load(sprintf('testknn1_data_%s.mat',options.metric));
end

fprintf('Test using %s metric\n',options.metric);
tic;
labelspred=knn1(data(1:nsamplestrain,:),labels(1:nsamplestrain,:),k,data(nsamplestrain+1:nsamples,:));
toc

if not(all(labelspred==labelspred_official))
    error('There is a problem with knn1');
else
    disp('knn1 test succeeded');
end