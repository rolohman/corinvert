function [ct_est,c0_est]=est_ct(d,Gi,cpmin)
d_orig=d;
[~,nd]=size(Gi);
nd=nd+1;
alpha=100;
cidl   = tril(ones(nd),-1)==1;
jnk=nan(nd);
jnk(cidl)=d;
jumps=mean(diff(jnk)','omitnan');
jumps(1)=jumps(2);

ct_est=zeros(1,nd-1);

[s,sortid]=sort(jumps);
sortid=sortid(s<0);

%mymax=sum(di.*w)/sum(w);

for i=sortid
    w=exp(alpha*d);
    dw=d.*w;
    id  = and(isfinite(d),Gi(:,i)==1); %box
    idn = and(isfinite(d),Gi(:,i)==0); %not box
    n   = sum(id);
    if(n>1)
        
        c1    = sum(dw(id))/sum(w(id));
        c2    = sum(dw(idn))/sum(w(idn));
       
    elseif(sum(id)==0)
        c1    = 0;
        c2    = 0;
    else
        c1    = max(d(id));
        c2    = min(d(id));
    end
    
    ct_est(i)=min(0,c1-c2);
    ct_est(i)=max(cpmin(i),ct_est(i));
    
    syn=Gi(:,i)*ct_est(i);
    d=d-syn;
end



