function mymax=mymax_mat(data,alpha)
%size data=[nsamp,nx], calcmax in first dim.
%i=isfinite(data);
%data(~i)=0;
mymax=sum(data.*exp(alpha*data),1,'omitnan')./sum(exp(alpha*data),1,'omitnan');
