function [gooddat,Grstat,Gistat]=find_bad(d,diags,permbad,Gr0,Gi0) 

%use these data
gooddat=~permbad;
gooddat(diags)=1;
gooddat(isnan(d))=0;

%invert for these model params
Grstat=sum(abs(Gr0(gooddat,:)),1)>2;
Gistat=sum(Gi0(gooddat,:),1)>1;
