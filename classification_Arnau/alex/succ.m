function ok = succ(yStd, clasfRow)
% Usage RateMatrix = succ(y, clasf)
%		y =  standard type class matrix
%		clasf =  vector with classnumber
%     RateMatrix return hits per class/class


[vecs,class]=size(yStd);
%Initialise matrix
% if class==2
%     1
% end
ok=zeros(class);

for ivec = 1:vecs
    %Current class  
    clv = find(yStd(ivec,:)==1);
    clf = clasfRow(ivec);
    
    ok(clv,clf) = ok(clv,clf) +1;
    
end


