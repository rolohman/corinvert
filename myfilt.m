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
xsum  = sum(windx);
ysum  = sum(windy);
windx = windx/xsum;
windy = windy/ysum;
ry    = floor(length(windy)/2);

z       = zeros(1,newnx);
fidi    = fopen(infile,'r');
fidm    = fopen(maskfile,'r');
fido    = fopen(outfile,'w');
fidc    = fopen(countfile,'w');

mask=fread(fidm,[newnx,newny],'integer*1');
mask=mask==1;

switch ftype
    case 1 %r4 phs
        in=fread(fidi,[newnx,newny],'real*4');
        in=exp(im*in);       
    case 2 %c8
        in=fread(fidi,[newnx*2, newny],'real*4');
        in=in(1:2:end,:)+im*in(2:2:end,:);   
    case 3 %unw r4
        in=fread(fidi,[newnx,newny],'real*4');
end

in(isnan(in))=0;
in(~mask)=0;

good= in~=0;

goodsum = conv2(windx,windy,good,'same');
datasum = conv2(windx,windy,in,'same');
out=datasum./goodsum;
%write count
fwrite(fidc,goodsum,'real*4');
switch(outtype)
    case 1 %just output filtered phase at all points
        fwrite(fido,angle(out),'real*4'); %1000pixel filtered product, at all pixels, even masked ones.
    case 2 %output filtered phase only at unmasked points, original elsewhere

    case 3 %output orig in - filtered phase %need to fix this so that it works for unmasked

    case 4 %output filtered, unwrapped input and output
        fwrite(fido,out,'real*4');
    case 5 %output unwrapped-filtered %need to fix this so that it works for unmasked

end


