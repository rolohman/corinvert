function mymax=mymax(data,alpha)
i=isfinite(data);
mymax=sum(data(i).*exp(alpha*data(i)))/sum(exp(alpha*data(i)));
