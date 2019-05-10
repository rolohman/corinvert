function mask_ztopo(topofile,file,type,nx,ny)

fidm=fopen(topofile,'r');
a=fread(fidm,[nx ny],'real*8');
switch type
    case 1 %r4
        fid1=fopen(file,'r');
        b=fread(fid1,[nx ny],'real*4');
        fclose(fid1);
        
        b(a<=0)=NaN;
        fid1=fopen(file,'w');
        fwrite(fid1,b,'real*4');
        fclose(fid1);
end