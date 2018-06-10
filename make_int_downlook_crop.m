function mycor_slcs_downlook(slcfile1,slcfile2,outfile,nx,ny,rx,ry,bx,by)

% outfile=cor file
% nx = width
% ny=length
% rx,ry half width of averaging window
% windowtype - 0=triangular weight, 1 = constant weight
% ampflag - 1: use amplitude file, 0: just phase

newnx=floor(nx/rx);
newny=floor(ny/ry);
nx2=bx(2)-bx(1)+1;
ny2=by(2)-by(1)+1;

rangevec=[0:newnx-1]*rx+1;
azvec=[0:newny-1]*ry+1;

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



        windx=ones(1,rx*2+1);
        windy=ones(1,ry*2+1);
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
       cpx=cpx1.*conj(cpx2);
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
    rea=circshift(rea,1);
    ima=circshift(ima,1);
    [a,count1]=fread(fid1,nx*2,'real*4');
    [b,count1]=fread(fid2,nx*2,'real*4');
    if(count1==nx*2)
        cpx1=a(1:2:end)+im*a(2:2:end);
        cpx2=b(1:2:end)+im*b(2:2:end);
        cpx=cpx1.*conj(cpx2);
        rea(1,:)=real(cpx);
        ima(1,:)=imag(cpx);
    else
        rea(1,:)=z;
        ima(1,:)=z;
    end
    
    
    mag    = sqrt(rea.^2+ima.^2);
    r2=rea;
    i2=ima;
    good=mag>0;
    r2(good)=r2(good)./mag(good);
    i2(good)=i2(good)./mag(good);
    
    r2     = windy*r2;
    i2     = windy*i2;
    
    rsum   = conv(r2,windx,'same');
    isum   = conv(i2,windx,'same');
    sm     = sqrt(rsum.^2+isum.^2);%length is one if completely coherent
    
    sm(isnan(sm))=0;
    if(ismember(j,azvec))
        fwrite(fid3,sm(rangevec),'real*4');
    end
end


fclose(fid3);
fclose(fid2);
fclose(fid1);
%system('cp amplitude.amp.rsc topophase_flat.cor.rsc')
