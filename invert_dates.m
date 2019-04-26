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

maskfile='cordir_VV/average.cor';
maskthresh=0.7;

tsdir='tsdir/';
if(~exist(tsdir,'dir'))
    mkdir(tsdir);
end
ddir=['tsdir/mskshift/']
if(~exist(ddir,'dir'))
    mkdir(ddir)
end

fid=fopen(maskfile,'r');
mask=fread(fid,[newnx newny],'real*4');
mask=mask>maskthresh;
fclose(fid);

for i=1:nd-1
    j=i+1;
    intdir=(['intdir' pol '/' dates(i).name '/']);
    ints(i).name=[dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
    ints(i).unw=[intdir ints(i).name];
    ints(i).msk=[ddir ints(i).name];
    
    fid=fopen(ints(i).unw,'r');
    tmp=fread(fid,[newnx,newny],'real*4');
    fclose(fid);
    tmp(~mask)=NaN;
    tmp=tmp-mean(tmp(mask),'omitnan');
    
    fid=fopen(ints(i).msk,'r');
    fwrite(fid,tmp,'real*4');
    fclose(fid);
    return
end
    
