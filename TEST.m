% acc = [];
% nMax = 100;
% for k = 1:nMax
%     KNN_k = fitcknn(train,gp,'NumNeighbors',k);
%     labels_k = predict(KNN_k,test2);
%     acc_k = sum(gp_test == labels_k)/size(labels_k,1);
%     acc = cat(1,acc,acc_k);
% end
% plot(1:nMax,acc);

acc = [];
nMax = 1000;
for k = 1:nMax
    idx = find(ftRank <= k);
    train_k = train(:,idx);
    test_k = test2(:,idx);
    svm_model = svmtrain(gp, train_k, '-t 0 -c 1');
    [svm_class, acc_k, dec] = svmpredict(gp_test, test2, svm_model);
    acc = cat(1,acc,acc_k(1));
end
plot(1:nMax,acc);
