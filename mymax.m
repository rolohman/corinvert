function mymax=mymax(data,alpha)
i=isfinite(data);
di=data(i);
w=exp(alpha*di);
mymax=sum(di.*w)/sum(w);

%Applications of lp-norms and their smooth approximations for gradient
%based learning vector quantization,Lange et al., 2014