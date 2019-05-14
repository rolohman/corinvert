function infill(int,msk,outint,filters,newnx,newny)
%filter multiple times and fill masked area
im=sqrt(-1);
for k=1:length(filters)
    filtname  = ['filt_' num2str(filters(k)) '.int'];
    filtwname = ['filt_' num2str(filters(k)) '.wgt'];
    myfilt(int,msk,filtname,filters(k),filters(k),newnx,newny,2,1,1,filtwname);
end

fid0 = fopen(msk,'r');
fid1 = fopen(int,'r');
for k=1:length(filters)
    filtname  = ['filt_' num2str(filters(k)) '.int'];
    filtwname = ['filt_' num2str(filters(k)) '.wgt'];
    
    fidf(k)=fopen(filtname,'r');
    fidw(k)=fopen(filtwname,'r');
end
fido = fopen(outint,'w');
outc = zeros(1,newnx*2);
for i=1:newny
    m  = fread(fid0,newnx,'integer*1');
    a  = fread(fid1,newnx,'real*4');
    
    wts=zeros(length(filters),newnx);
    fis=wts;
    for j=1:length(filters)
        fis(j,:)=fread(fidf(k),newnx,'real*4');
        wts(j,:)=fread(fidw(k),newnx,'real*4');
    end
       
    a   = exp(im*a);
    fis = exp(im*fis);
             
    m1 = m==1;
            
    %a, unless masked, then combine filtered products, tapered.
    out(m1)           = a(m1);
 
    tmpwgt=[wgts(1,:);diff(wts,1,1)]; %wedges
    tot=sum(tmpwgt,1);
    tmpwgt=tmpwgt./repmat(tot,length(filters,1));
           
    sums=sum(tmpwgt.*fis,1,'omitnan');
    out(~m1) = sums(~m1);
    out(isnan(a)) = NaN;

    outc(1:2:end)     = cos(out);
    outc(2:2:end)     = sin(out);
    outc(isnan(outc)) = 0;
    fwrite(fido,outc,'real*4');
end
        
 