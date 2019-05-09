function mask_ztopo(topofile,file,type,nx,ny)

fidm=fopen(topofile,'r');
fid1=fopen(file,'r');


a=fread(fidm,[nx ny],'real*8');
switch type
    case 1 %r4
        b=fread(fid1,[nx ny],'real*4');
        b(a<=0)=0;
        fwrite(fido1,b,'real*4');
end