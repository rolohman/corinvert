function [ct_est,c0_est]=est_ct_c0(d,Gi,cpmin)
d_orig=d;
[ni,nd]=size(Gi);
nd=nd+1;
alpha=100;
cidl   = tril(ones(nd),-1)==1;
jnk=nan(nd);
jnk(cidl)=d;
jumps=median(diff(jnk)','omitnan');

ct_est=zeros(1,nd-1);
c0_est=mymax(d,alpha);

notdone=1;

[s,sortid]=sort(jumps);
sortid=sortid(s<0.01);

for i=sortid

        id  = and(isfinite(d),Gi(:,i)==1); %box
        idn=~id;
        n   = sum(id);
        if(n>1)
            c1    = mymax(d(id),alpha);
            c2    = mymax(d(idn),alpha);
        elseif(sum(id)==0)
            c1    = 0;
            c2    = 0;
        else
            c1    = max(d(id));
            c2    = min(d(id));
        end
        
        ct_est(i)=min(0,c1-c2);
        ct_est(i)=max(cpmin(i),ct_est(i));
        
        syn=Gi*ct_est';
        d=d_orig-syn;
    
end


alpha=25;
c0_est=mymax(d,alpha);


