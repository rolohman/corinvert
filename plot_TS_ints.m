function [unw, msk,hp,hpfx,dem]=plot_TS_ints(x,y)
pol  = '_VV';
ddir = ['dates' pol '/'];
decide_ints_stack

for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).dn      = datenum(dates(i).name,'YYYYmmdd');
 end
dn=[dates.dn];
dn=dn-min(dn);

clear ints
for i=1:nd-1
    j=i+1;
    ints(i).dir   = (['intdir' pol '/' dates(i).name '/']);
    ints(i).name  = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).unw   = [ints(i).dir ints(i).name '.unw'];
    ints(i).msk   = [ints(i).dir ints(i).name '.msk'];
    ints(i).hp    = [ints(i).dir ints(i).name '_highpass.unw'];
    ints(i).hpfx  = [ints(i).dir ints(i).name '_highpass_fix.unw'];
end   



fid=fopen(['demerr.r4'],'r');
fseek(fid,((y-1)*newnx+x-1)*4,-1);
dem=fread(fid,1,'real*4');


for i=1:nd-1
    fid=fopen(ints(i).msk,'r');
    fseek(fid,((y-1)*newnx+x-1),-1);
    msk(i)=fread(fid,1,'integer*1');
    fid=fopen(ints(i).unw,'r');
    fseek(fid,((y-1)*newnx+x-1),-1);
    unw(i)=fread(fid,1,'real*4');
    fid=fopen(ints(i).hp,'r');
    fseek(fid,((y-1)*newnx+x-1),-1);
    hp(i)=fread(fid,1,'real*4');
    fid=fopen(ints(i).hpfx,'r');
    fseek(fid,((y-1)*newnx+x-1),-1);
    hpfx(i)=fread(fid,1,'real*4');
end

