function [crm,cp]=flatten(cmed,alpha,cp,cpmin)
nd=length(cmed);
for i=1:nd
    crm1(i)=mymax(cmed(i:end),alpha);
    crm2(i)=mymax(cmed(1:i),alpha);
end

crm1=cummax(crm1,'reverse');
crm2=cummax(crm2);
crm=crm1+crm2;
crm(end)=crm(end-1);
crm(1)=crm(2);

cdif=diff(crm);

cdif=-abs(cdif);

if(nargin>2)
    testcp=cp+cdif;
    badi=testcp<cpmin;
    if(sum(badi))
        
        badi=find(badi);
        bdif=cpmin(badi)-testcp(badi);
        cdif(badi)=cdif(badi)+bdif;
        
        for i=1:length(badi)
            crm(badi(i)+1:end)=crm(badi(i)+1:end)+bdif(i);
        end
    end
    
    
    cp=cp+cdif;
end