function [res,mags,shift,G]=expfun(model,x,y)
% global timesall
% timesall(end+1,:)=model;
[nd,nr]=size(x);
 
ts=repmat(model,nd,1);

%options=optimset('display','none');
G=exp(-x./ts);
G(x==0)=0;
G(:,end+1)=ones(nd,1);

[newmod,rsn,res]=lsqlin(G,y);
%[newmod,rsn,res]=lsqlin(G,y,[],[],[],[],[zeros(1,nr) -1],max(y)*2*ones(1,nr+1),[],options);
mags=newmod(1:end-1);
shift=newmod(end);
