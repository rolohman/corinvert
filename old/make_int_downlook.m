function make_int_downlook(slcfile1,slcfile2,outfile,nx,ny,rx,ry,windowtype)

% outfile=cor file
% nx = width
% ny=length
% rx,ry half width of averaging window
% windowtype - 0=triangular weight, 1 = constant weight
% ampflag - 1: use amplitude file, 0: just phase

newnx=floor(nx/rx);
newny=floor(ny/ry);
rangevec=[1:newnx]*rx-floor(rx/2);
azvec=[1:newny]*ry-floor(ry/2);
if(exist(slcfile1,'file'))
    fid1=fopen(slcfile1,'r');
else
    disp([slcfile1 ' does not exist'])
    return
end
if(exist(slcfile2,'file'))
    fid2=fopen(slcfile2,'r');
else
    disp([slcfile2 ' does not exist'])
    return
end

fid3=fopen(outfile,'w');
im=sqrt(-1);

switch windowtype
    case 0
        windx=[1/(rx+1):1/(rx+1):1 (1-1/(rx+1)):-1/(rx+1):1/(rx+1)];
        windy=[1/(ry+1):1/(ry+1):1 (1-1/(ry+1)):-1/(ry+1):1/(ry+1)];
    case 1
        windx=ones(1,rx*2+1);
        windy=ones(1,ry*2+1);
    case 2
        windx=exp(-(-rx:rx).^2/2/(rx/2)^2);
        windy=exp(-(-ry:ry).^2/2/(ry/2)^2);
end
windx=windx/sum(windx);
windy=windy/sum(windy);
ry=floor(length(windy)/2);

z   = zeros(1,nx);
rea = zeros(ry*2+1,nx);
ima = zeros(ry*2+1,nx);


%load first ry lines
for j=1:ry%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rea=circshift(rea,1);
    ima=circshift(ima,1);
    
    [a,count1]=fread(fid1,nx*2,'real*4');
    [b,count1]=fread(fid2,nx*2,'real*4');
     if(count1==nx*2)
        cpx1=a(1:2:end)+im*a(2:2:end);
        cpx2=b(1:2:end)+im*b(2:2:end);
          
        a1=abs(cpx1);
        a2=abs(cpx1);
        cpx=cpx1.*conj(cpx2);
        am=sqrt(a1.*a2);
        goodid=am>0;
        cpx(goodid)=cpx(goodid)./am(goodid);
        rea(1,:)=real(cpx);
        ima(1,:)=imag(cpx);
    else
        rea(1,:)=z;
        ima(1,:)=z;

    end    
end

%now go to end (passing end by ry lines, filling in with zeros. "active"
%line is at ry+1th row

for j=1:ny
    cpxvec=zeros(1,2*newnx);
    rea=circshift(rea,1);
    ima=circshift(ima,1);
    [a,count1]=fread(fid1,nx*2,'real*4');
    [b,count1]=fread(fid2,nx*2,'real*4');

    if(count1==nx*2)
        cpx1=a(1:2:end)+im*a(2:2:end);
        cpx2=b(1:2:end)+im*b(2:2:end);
        
        a1=abs(cpx1);
        a2=abs(cpx2);
        cpx=cpx1.*conj(cpx2);
        am=sqrt(a1.*a2);
        goodid=am>0;
        cpx(goodid)=cpx(goodid)./am(goodid);
        
        rea(1,:)=real(cpx);
        ima(1,:)=imag(cpx);
    else
        rea(1,:)=z;
        ima(1,:)=z;
        am=z;
    end
    
    mag  = sqrt(rea.^2+ima.^2);
    
    r2   = rea;
    i2   = ima;
    r2   = windy*r2;
    i2   = windy*i2;
    m2   = windy*mag;
    
    r2(isnan(r2))=0;
    i2(isnan(i2))=0;
    m2(isnan(m2))=1;
    rsum = conv(r2,windx,'same');
    isum = conv(i2,windx,'same');
    msum = conv(m2,windx,'same');
    cpx3 = rsum+im*isum;
    cpx3 = cpx3./msum;
    sm   = abs(cpx3);
    sm(isnan(sm))=0;
    cpxvec(1:2:end)=rsum(rangevec);
    cpxvec(2:2:end)=isum(rangevec);
    if(ismember(j,azvec))
        fwrite(fid3,cpxvec,'real*4');
    end
end


fclose(fid3);
fclose(fid2);
fclose(fid1);

