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
    fid2=fopen(dates(i).fixunw,'r');
    fseek(fid2,((y-1)*newnx+x-1)*4,-1);
    def2(i)=fread(fid,1,'real*4');
end

synth=intercept+slope*dn;
figure
plot(dn,def,'.-')
hold on
plot(dn,def2,'.-')
plot(dn,synth,'r')
