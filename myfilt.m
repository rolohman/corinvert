function myfilt(infile,maskfile,outfile,rx,ry,newnx,newny,ftype,outtype)
im=sqrt(-1);
windx=ones(1,rx*2+1);
windy=ones(1,ry*2+1);
windx=windx/sum(windx);
windy=windy/sum(windy);


z       = zeros(1,newnx);
fidi    = fopen(infile,'r');
fidm    = fopen(maskfile,'r');
fido    = fopen(outfile,'w');

mask=fread(fidm,[newnx,ry],'integer*1');
mask=mask==1;
mask=[flipud(mask');false(ry+1,newnx)];


switch ftype
    case 1 %r4
        in=fread(fidi,[newnx,ry],'real*4');
        in = exp(im*in);
    case 2 %int
        in=fread(fidi,[newnx*2, ry],'real*4');
        in=in(1:2:end,:)+im*in(2:2:end,:);
end

in = [flipud(in');nan(ry+1,newnx)];

in(isnan(in))=0;
in(~mask)=0;

%now go to end (passing end by ry lines, filling in with zeros. "active"
%line is at ry+1th row

for j=1:newny
    mask=circshift(mask,1);
    [a,count1]=fread(fidm,newnx,'integer*1');
    a(isnan(a))=false;
    if(count1==newnx)
        
        mask(1,:)=a;
    else
        mask(1,:)=z;
    end
    
    
    in=circshift(in,1);
    switch ftype
        case 1 %r4
            a=fread(fidi,[newnx,ry],'real*4');
            a = exp(im*a);
        case 2 %int
            a=fread(fidi,[newnx*2, ry],'real*4');
            a=a(1:2:end,:)+im*a(2:2:end,:);
    end
    a(isnan(a))=0;
    if(count1==newnx)
        a(~mask(1,:))=0;
        in(1,:)=a;
    else
        in(1,:)=z;
    end
    good = in~=0;
    a    = good;
    c    = in;
    
    a    = windy*a;
    c    = windy*c;
    
    asum = conv(a,windx,'same');
    csum = conv(c,windx,'same');
    out  = csum./asum;
    
    switch(outtype)
        case 1 %just output filtered phase at all points
            fwrite(fido,angle(out),'real*4'); %1000pixel filtered product, at all pixels, even masked ones.
        case 2 %output filtered phase only at unmasked points
            orig=in(:,ry+1);
            msk = mask(:,ry+1);
            out(msk)=orig(msk);
            fwrite(fido,angle(out),'real*4');
        case 3 %output orig in - filtered phase
            orig=in(:,ry+1);
            dif=orig*conj(out);
            fwrite(fido,angle(dif),'real*4');
    end
    
end
