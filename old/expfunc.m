function [res,synth,mags]=expfunc(mod,dn,dnr,data)
ndates=length(dnr);
%de=max(0.03,0.25*data);
%de=1./de;

Gs=zeros(length(dn),ndates);
for i=1:ndates
    Gs(:,i)=exp(-(dn-dnr(i))/mod(i));
    Gs(dn<dnr(i),i)=0;
end

   
G=[Gs];


mags=lsqnonneg(G,data);

synth=G*mags;

res=data-synth;
%res=de.*res;


