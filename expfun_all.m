function [res]=expfun_all(model,x,y)
% global timesall
% timesall(end+1,:)=model;
[nd,nr]=size(x);

ts=repmat(model(1:nr),nd,1);
%options=optimset('display','none');
G=exp(-x./ts);
G(x==0)=0;
G(:,end+1)=ones(nd,1);

fwd=G*model(nr+1:end)';
res=y-fwd;


