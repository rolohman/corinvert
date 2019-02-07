function mymax=mymax(data,alpha)
i=isfinite(data);
di=data(i);
w=exp(alpha*di);
mymax=sum(di.*w)/sum(w);
