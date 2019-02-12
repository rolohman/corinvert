decide_ints_stack


difdir='merged/diffint';
if(~exist(difdir,'dir'))
    mkdir(difdir)
end

for i=1:nd-1
    slc1=[dates(i).name  '/' dates(i).name '.slc.full'];
    slc2=[dates(i+1).name '/' dates(i+1).name '.slc.full'];
    out=[difdir '/' dates(i).name '-' dates(i+1).name 'diff.int'];
    if(~exist(out,'file'))
        command=['imageMath.py -e=''a*conj(b)*conj(c)*f'' -t cfloat -o tmp --a=merged/SLC_VV/' slc1 ' --b=merged/SLC_VV/' slc2 ' --c=merged/SLC_VH/' slc1 ' --f=merged/SLC_VH/' slc2];
        system(command);
       command=['imageMath.py -e=''a/abs(a)*sqrt(sqrt(abs(a)))'' -t cfloat -o ' out ' --a=tmp']; 
           system(command);
    else
        disp([out ' already made'])
    end
end

%now look at avg.
for i=1:nd-1
    out=[difdir '/' dates(i).name '-' dates(i+1).name 'diff.int'];
    fid(i)=fopen(out,'r');
end
im=sqrt(-1);
fo=fopen('meanphs2.int','w');
for i=1:ny
    for j=1:nd-1
        tmp=fread(fid(j),nx*2,'real*4');
        cpx(j,:)=tmp(1:2:end)+im*tmp(2:2:end);
    end
    cpx=cpx./abs(cpx); 
    cpxm=mean(cpx,1,'omitnan');
    %cpxm=mean(cpx,1,'omitnan')./mean(abs(cpx),1,'omitnan');
    output=zeros(1,nx*2);
    output(1:2:end)=real(cpxm);
    output(2:2:end)=imag(cpxm);
fwrite(fo,output,'real*4');
end
fclose('all');
