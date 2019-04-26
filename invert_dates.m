pol='_VV';
decide_ints_stack
slcdir=['merged/SLC' pol '/'];
skips=1;  %1=sequential, larger=longer time pairs
for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
    dates(i).dn      = datenum(dates(i).name,'YYYYmmdd');
end
dn=[dates.dn];

maskfile='wgtdir/average.1alks_3rlks.amp';
maskthresh=0.14;

tsdir='tsdir/';
if(~exist(tsdir,'dir'))
    mkdir(tsdir);
end
ddir=['tsdir/mskshift/'];
if(~exist(ddir,'dir'))
    mkdir(ddir)
end
fdir=['tsdir/filt/'];
if(~exist(fdir,'dir'))
    mkdir(fdir)
end

fid=fopen(maskfile,'r');
mask=fread(fid,[newnx newny],'real*4');
mask=mask>maskthresh;
fclose(fid);
geomfile='merged/geom_master/hgt.rdr.1alks_3rlks.full';
fid=fopen(geomfile,'r');
for i=1:newny
    a=fread(fid,newnx,'real*8');
    mask(a<0,i)=0;
end
fclose(fid);

for i=1:nd-1
    j=i+1;
    intdir=(['intdir' pol '/' dates(i).name '/']);
    ints(i).name=[dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
    ints(i).unw=[intdir ints(i).name];
    ints(i).msk=[ddir ints(i).name];
    ints(i).filt=[fdir ints(i).name];
    
end

for i=1:nd-1
    if(~exist(ints(i).msk,'file'))
        fid=fopen(ints(i).unw,'r');
        tmp=fread(fid,[newnx,newny],'real*4');
        fclose(fid);
        
        tmp(~mask)=NaN;
        tmp=tmp-mean(tmp(mask),'omitnan');
        fid=fopen(ints(i).msk,'w');
        fwrite(fid,tmp,'real*4');
        fclose(fid);
    end
    return
end

ry=500;
windx=ones(1,ry*2+1);
windy=ones(1,ry*2+1);
windx=windx/sum(windx);
windy=windy/sum(windy);
for i=1:nd-1
    if(~exist(ints(i).filt,'file'))
        
        z       = zeros(1,newnx);
      
        fidi    = fopen(ints(i).msk,'r');
        fido    = fopen(ints(i).filt,'w');
        
        mskunw=fread(fidi,[newnx,ry],'real*4');
        mskunw = [flipud(mskunw');nan(ry+1,newnx)];
        mskunw(isnan(mskunw))=0;
        
        %now go to end (passing end by ry lines, filling in with zeros. "active"
        %line is at ry+1th row
        
        for j=1:newny
            
            mskunw=circshift(mskunw,1);
            [a,count1]=fread(fidi,newnx,'real*4');
            a(isnan(a))=0;
            if(count1==newnx)
                mskunw(1,:)=a;
            else
                mskunw(1,:)=z;
            end
            good=mskunw~=0;
            a   = good;
            c   = mskunw;
            
            a   = windy*a;
            c   = windy*c;
            
            asum = conv(a,windx,'same');
            csum = conv(c,windx,'same');
            out = csum./asum;
            
            fwrite(fido,out,'real*4'); %1000pixel filtered product, at all pixels, even masked ones.
            
        end
        fclose(fidi);
        fclose(fido);
    end
 
end
