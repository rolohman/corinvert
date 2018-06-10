function [ct_est,c0_est,synth,cres]=est_ct_c0_smoothmax(data,id1,id2,Gi)
[ni,nd]=size(Gi);
nd=nd+1;

alpha=100;
c0_est=mymax(data,alpha);

data=data/c0_est;
c1=ones(1,nd);

for i=1:nd-1
    id  = and(isfinite(data),and(id1<=i,id2>i)); %box
    idb = and(isfinite(data),id2==i); %column before box
    idf = and(isfinite(data),id2==i+1); %first column in box
    n   = sum(id);
    nb  = sum(idb);
    nf  = sum(idf);
    if(n>5)
        c1(i+1)    = mymax(data(id),alpha);
    elseif(sum(id)==0)
        c1(i+1)    = c1(i);
    else
        c1(i+1)    = max(data(id));
    end
    
    if(and(nb>5,nf>5))
        cb         = mymax(data(idb),alpha);
        cf         = mymax(data(idf),alpha);
        cvarb      = ((0.3083*(1-cb)-.2*(1-cb).^2))^2/nb;
        cvarf      = ((0.3083*(1-cf)-.2*(1-cf).^2))^2/nf;
        cdifstd(i) = sqrt(cvarb+cvarf);
        cdif(i)    = cb-cf;
    elseif(or(nb==0,nf==0))
        cdifstd(i) = inf;
        cdif(i)    = 0;
    else
        cb         = median(data(idb));
        cf         = median(data(idf));
        cvarb      = ((0.3083*(1-cb)-.2*(1-cb).^2))^2/nb;
        cvarf      = ((0.3083*(1-cf)-.2*(1-cf).^2))^2/nf;
        cdifstd(i) = sqrt(cvarb+cvarf);
        cdif(i)    = cb-cf;
    end
    
end
c1(1)=c1(2);

ct_est = 1+diff(c1);
ct_est(or(ct_est>1,cdif<cdifstd))=1;

synth=exp(Gi*log(ct_est'));
cres=data./synth/c0_est;

