function [crm,cp]=flatten_back(cmed,alpha,cp,cpmin)
nd=length(cmed);
for i=1:nd
    crm(i)=mymax(cmed(i:end),alpha);
end

crm=cummax(crm,'reverse');

crm(end)=crm(end-1);
crm(1)=crm(2);
cdif=-abs(diff(crm));
cdif(isnan(cdif))=0;

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