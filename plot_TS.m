function [def, slope, intercept, res,def2]=plot_TS(x,y)
pol  = '_VV';
ddir = ['dates' pol '/'];
decide_ints_stack

for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).dn      = datenum(dates(i).name,'YYYYmmdd');
    dates(i).unw     = [ddir dates(i).name '_simple.unw'];
    dates(i).fixunw     = [ddir dates(i).name '_simple_fix.unw'];
    dates(i).infill  = [ddir dates(i).name '_infill.unw'];
end
dn=[dates.dn];
dn=dn-min(dn);

clear ints
for i=1:nd-1
    j=i+1;
    ints(i).dir   = (['intdir' pol '/' dates(i).name '/']);
    ints(i).name  = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).int   = [ints(i).dir ints(i).name '.int'];
    ints(i).msk   = [ints(i).dir ints(i).name '.msk'];
end
   



fid=fopen([ddir 'simpleslope.r4'],'r');
fseek(fid,((y-1)*newnx+x-1)*4,-1);
slope=fread(fid,1,'real*4');

fid=fopen([ddir 'simpleint.r4'],'r');
fseek(fid,((y-1)*newnx+x-1)*4,-1);
intercept=fread(fid,1,'real*4');

fid=fopen([ddir 'simpleres.r4'],'r');
fseek(fid,((y-1)*newnx+x-1)*4,-1);
res=fread(fid,1,'real*4');


for i=1:nd
    fid=fopen(dates(i).unw,'r');
    fseek(fid,((y-1)*newnx+x-1)*4,-1);
    def(i)=fread(fid,1,'real*4');
    fid=fopen(dates(i).fixunw,'r');
    fseek(fid,((y-1)*newnx+x-1)*4,-1);
    def2(i)=fread(fid,1,'real*4');
end

for i=1:nd-1
    fid=fopen(ints(i).msk,'r');
    fseek(fid,((y-1)*newnx+x-1),-1);
    msk(i)=fread(fid,1,'integer*1');
end

synth=intercept+slope*dn;
figure
plot(dn,def,'.-')
hold on
plot(dn,def2,'.-')
plot(dn,synth,'r')
for i=1:nd
    if(~msk(i))
        plot(dn(i),def(i),'ko')
    end
end