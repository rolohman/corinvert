function make_intcor_anyposting(slcfile1,slcfile2,corfile,intfile,nx,ny,px,py,rx,ry,windowtype,wgtfile)

% outfile=cor file
% nx = width
% ny=length
% rx,ry half width of averaging window
% windowtype - 0=triangular weight, 1 = constant weight
% ampflag - 1: use amplitude file, 0: just phase
% wgtfile: must be 0-1.

newnx=floor(nx/px);
newny=floor(ny/py);
rangevec=[1:newnx]*px-floor(px/2);
azvec=[1:newny]*py-floor(py/2);
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
if(exist('wgtfile','var'))
    if(wgtfile)
        uwgt=1;
        if(exist(wgtfile,'file'))
            fidw=fopen(wgtfile,'r');
        else
            disp([wgtfile ' does not exist']);
            uwgt=0;
        end
    else
        uwgt=0;
        disp('not using external wgt file')
    end
else
    uwgt=0;
    disp('not using external wgt file')
end

fid3=fopen(corfile,'w');
fid4=fopen(intfile,'w');
im=sqrt(-1);

switch windowtype
    case 0
        windx=[1/(rx+1):1/(rx+1):1 (1-1/(rx+1)):-1/(rx+1):1/(rx+1)];
        windy=[1/(ry+1):1/(ry+1):1 (1-1/(ry+1)):-1/(ry+1):1/(ry+1)];
    case 1
        windx=zeros(1,rx*2+1); windx((1:rx)+ceil(rx/2))=1;
        windy=zeros(1,ry*2+1); windy((1:ry)+ceil(ry/2))=1;
    case 2
        windx=exp(-(-rx*2:rx*2).^2/2/(rx/2)^2);
        windy=exp(-(-ry*2:ry*2).^2/2/(ry/2)^2);
end
windx=windx/sum(windx);
windy=windy/sum(windy);
ry=floor(length(windy)/2);

z    = zeros(1,nx);
slc1  = zeros(ry*2+1,nx);
slc2  = slc1;
wgts  = slc1;

%load first ry lines
for j=1:ry%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    slc1=circshift(slc1,1);
    slc2=circshift(slc2,1);
    [a,count1]=fread(fid1,nx*2,'real*4');
    [b,count1]=fread(fid2,nx*2,'real*4');
    if(count1==nx*2)
        cpx1=a(1:2:end)+im*a(2:2:end);
        cpx2=b(1:2:end)+im*b(2:2:end);
        
        slc1(1,:)=cpx1;
        slc2(1,:)=cpx2;
    else
        slc1(1,:)=z;
        slc2(1,:)=z;
    end
    if(uwgt)
        wgts=circshift(wgts,1);
        [w,count]=fread(fidw,nx,'real*4');
        w(~isfinite(w))=0;
        if(count==nx)
            wgts(1,:)=w;
        else
            wgts(1,:)=z;
        end
    end
end
%now go to end (passing end by ry lines, filling in with zeros. "active"
%line is at ry+1th row

for j=1:ny
    
    slc1=circshift(slc1,1);
    slc2=circshift(slc2,1);
    [a,count1]=fread(fid1,nx*2,'real*4');
    [b,count1]=fread(fid2,nx*2,'real*4');
    if(count1==nx*2)
        cpx1=a(1:2:end)+im*a(2:2:end);
        cpx2=b(1:2:end)+im*b(2:2:end);
        
        slc1(1,:)=cpx1;
        slc2(1,:)=cpx2;
    else
        slc1(1,:)=z;
        slc2(1,:)=z;
    end
    if(uwgt)
        wgts=circshift(wgts,1);
        [w,count]=fread(fidw,nx,'real*4');
        w(~isfinite(w))=0;
        if(count==nx)
            wgts(1,:)=w;
        else
            wgts(1,:)=z;
        end
        %wgts(wgts<0.2)=0;
    end
    if(ismember(j,azvec))
        if(uwgt)
            a   = slc1.*conj(slc1).*wgts;
            b   = slc2.*conj(slc2).*wgts;
            c   = slc1.*conj(slc2).*wgts;
        else
            a   = slc1.*conj(slc1);
            b   = slc2.*conj(slc2);
            c   = slc1.*conj(slc2);
        end
        a   = windy*a;
        b   = windy*b;
        c   = windy*c;
        %disp([j length(isnan(a)) length(isnan(b)) length(isnan(c))])
        asum = conv(a,windx,'same');
        bsum = conv(b,windx,'same');
        csum = conv(c,windx,'same');
        cpx3 = csum./sqrt(asum.*bsum);
        cpx3 = cpx3(rangevec);
        
        fwrite(fid3,abs(cpx3),'real*4'); %cor
        fwrite(fid4,angle(cpx3),'real*4'); %int
    end    
end

fclose(fid4);
fclose(fid3);
fclose(fid2);
fclose(fid1);

