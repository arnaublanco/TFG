function ySt=yRow2yStd(yin)
% Usage ySt=yRow2yStd(yin)
%    This fucntion converts Rows labels style 
%    to ISPG Standard binary label style
%     
%  Alex Perera   ISPG/eSenseSystems 
%  15 May 2001

[vecs isone]= size(yin);

if isone~=1,
   error('yRow2yStd: not column format');
end;
ySt=[];
itd = itmark(0,vecs);
for k=1:vecs,
   ySt(k,yin(k))=1;
   itd = itmark(itd,vecs);
end
itd = itmark(-1,vecs);
   
