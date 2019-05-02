function myfilt(infile,maskfile,outfile,rx,ry,newnx,newny,windowtype,ftype,outtype,countfile)
im=sqrt(-1);



switch windowtype
    case 0
        windx=[1/(rx+1):1/(rx+1):1 (1-1/(rx+1)):-1/(rx+1):1/(rx+1)];
        windy=[1/(ry+1):1/(ry+1):1 (1-1/(ry+1)):-1/(ry+1):1/(ry+1)];
    case 1
        windx=zeros(1,rx*2+1); windx((1:rx)+ceil(rx/2))=1;
        windy=zeros(1,ry*2+1); windy((1:ry)+ceil(ry/2))=1;
    case 2 %gaussian
        windx=exp(-(-rx*2:rx*2).^2/2/(rx/2)^2);
        windy=exp(-(-ry*2:ry*2).^2/2/(ry/2)^2);
    case 3 %exp
        windx=exp(-abs(-rx*2:rx*2)/rx*3);
        windy=exp(-abs(-ry*2:ry*2)/ry*3);
end
xsum=sum(windx);
ysum=sum(windy);
windx=windx/xsum;
windy=windy/ysum;
ry=floor(length(windy)/2);

z       = zeros(1,newnx);
fidi    = fopen(infile,'r');
fidm    = fopen(maskfile,'r');
fido    = fopen(outfile,'w');
fidc    = fopen(countfile,'w');

mask=fread(fidm,[newnx,ry],'integer*1');
mask=mask==1;
mask=[flipud(mask');false(ry+1,newnx)];


switch ftype
    case 1 %r4 phs
        in=fread(fidi,[newnx,ry],'real*4');
        in=exp(im*in);       
    case 2 %c8
        in=fread(fidi,[newnx*2, ry],'real*4');
        in=in(:,1:2:end)+im*in(:,2:2:end);   
    case 3 %unw r4
        in=fread(fidi,[newnx,ry],'real*4');
end

in = [flipud(in.');nan(ry+1,newnx)];
in(isnan(in))=0;
in(~mask)=0;


%now go to end (passing end by ry lines, filling in with zeros. "active"
%line is at ry+1th row

for j=1:newny
    mask=circshift(mask,1);
    [a,count1]=fread(fidm,newnx,'integer*1');
    a=a==1;
    a(isnan(a))=false;
    if(count1==newnx)       
        mask(1,:)=a;
    else
        mask(1,:)=z;
    end
    
    
    in=circshift(in,1);
    switch ftype
        case 1 %r4
            a=fread(fidi,newnx,'real*4');
            a = exp(im*a);
        case 2 %int
            a=fread(fidi,newnx*2,'real*4');
            a=a(1:2:end)+im*a(2:2:end);
        case 3 %unw r4
            a=fread(fidi,newnx,'real*4');
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
    out  = csum.*asum;

    %write count
    fwrite(fidc,asum,'real*4');
    switch(outtype)
        case 1 %just output filtered phase at all points
            fwrite(fido,angle(out),'real*4'); %1000pixel filtered product, at all pixels, even masked ones.
        case 2 %output filtered phase only at unmasked points, original elsewhere
            orig     = in(ry+1,:);
            msk      = mask(ry+1,:);
            out(msk) = orig(msk);
            fwrite(fido,angle(out),'real*4');
        case 3 %output orig in - filtered phase %need to fix this so that it works for unmasked
            orig     = in(:,ry+1);
            dif      = orig*conj(out);
            fwrite(fido,angle(dif),'real*4');
        case 4 %output filtered, unwrapped input and output
            fwrite(fido,out,'real*4');
        case 5 %output unwrapped-filtered %need to fix this so that it works for unmasked
            orig     = in(ry+1,:);
            fwrite(fido,orig-out,'real*4');
    end
    
end
