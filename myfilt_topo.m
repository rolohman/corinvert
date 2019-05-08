function myfilt_topo(infile,avgwgtfile,corfile,outfile,rx,ry,newnx,newny,demerrfile)
im=sqrt(-1);

fid       = fopen(avgwgtfile,'r');
mask1     = fread(fid,[newnx newny],'real*4');
fid       = fopen(corfile,'r');
mask2     = fread(fid,[newnx newny],'real*4');
mask      = and(mask1>0.1,mask2>0.4);

%gaussian
windx=exp(-(-rx*2:rx*2).^2/2/(rx/2)^2);
windy=exp(-(-ry*2:ry*2).^2/2/(ry/2)^2);

xsum  = sum(windx);
ysum  = sum(windy);
windx = windx/xsum;
windy = windy/ysum;
ry    = floor(length(windy)/2);

z       = zeros(1,newnx);
fidi    = fopen(infile,'r');

fido    = fopen(outfile,'w');
fidc    = fopen(demerrfile,'w');

[bp,intbp]=get_baselines;

in=fread(fidi,[newnx,newny],'real*4');
in=exp(im*in);

orig=in;
in(isnan(in))=0;
in(~mask)=0;

good= in~=0;

goodsum = conv2(windx,windy,good,'same');
datasum = conv2(windx,windy,in,'same');
filtered=datasum./goodsum;

flattened=in.*conj(filtered);

fid=fopen(outfile);
fwrite(fid,angle(flattened),'real*4');
fclose('all');



