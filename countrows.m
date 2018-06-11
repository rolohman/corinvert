function [permbad,medrow,ratio]=countrows(id1,id2,cors,mncorl,mncor)
[ni,nx]=size(cors);
nd=max(id2);
medrow=zeros(ni,nx);
medrect=zeros(ni,nx);
nbef=zeros(ni,nx);
cbef=zeros(ni,nx);

d=log(cors);


for k=1:ni    
    recid=and(id1<=id1(k),id2>=id2(k));
    afids=and(id1==id1(k),id2>=id2(k));
    befids=and(id1==id1(k),id2<=id2(k));
    medrow(k,:)=median(d(afids,:),'omitnan');
    medrect(k,:)=median(d(recid,:),'omitnan');
     nbef(k,:)=sum(befids);
    cbef(k,:)=sum(d(befids,:)>mncorl,1,'omitnan');  
end
ratio=cbef./nbef;


permbad = and(medrow<log(mncor*2),ratio<0.85); 
permbad = and(permbad,medrect<mncorl);
permbad = or(permbad,isnan(cors));
%now set permbad=true if any previous point on the row was true
% return
% badids=find(permbad);
% for k=badids'
%     row=k+[id2(k):nd]-id2(k);
%     permbad(row)=1;
% end

%but always use diagonals
%permstat(id2==id1+1)=0;