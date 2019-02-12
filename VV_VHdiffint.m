decide_ints_stack
slcdir=['merged/SLC' pol '/'];
for i=1:nd
    dates(i).name    = files(i).name(1:8);    
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
end

difdir='merged/diffint';
if(~exist(difdir,'dir'))
    mkdir(difdir)
end

for i=1:nd-1
    slc1=[dates(i).name  '/' dates(i).name '.slc.full'];
    slc2=[dates(i+1).name '/' dates(i+1).name '.slc.full'];
    out=[difdir '/' dates(i).name '-' dates(i+1).name 'diff.int'];
    if(~exist(out,'file'))
        command=['imageMath.py -e=''a*conj(b)*conj(c)*f'' -t cfloat -o ' out ' --a=SLC_VV/' slc1 ' --b=SLC_VV/' slc2 ' --c=SLC_VH/' slc1 ' --f=SLC_VH/' slc2];
        system(command)
    else
        disp([out ' already made'])
    end
end

