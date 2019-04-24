%parpool(8)
pol='_VV';
decide_ints_stack
slcdir=['merged/SLC' pol '/'];
skips=1;  %1=sequential, larger=longer time pairs
for i=1:nd
    dates(i).name    = files(i).name(1:8);    
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
end

slcdirVH=['merged/SLC_VH/'];
for i=1:nd   
    dates(i).slcVH     = [slcdirVH dates(i).name '/' dates(i).name '.slc.full'];
end

if(~exist(['cordir' pol]))
    mkdir(['cordir' pol])
end
if(~exist(['intdir' pol]))
    mkdir(['intdir' pol])
end
if(~exist(['wgtdir']))
    mkdir(['wgtdir'])
end
for i=1:nd-1
    cordir=(['cordir' pol '/' dates(i).name '/']);
    intdir=(['intdir' pol '/' dates(i).name '/']);
    if(~exist(cordir,'dir'))
        mkdir(cordir)
    end
    if(~exist(intdir,'dir'))
        mkdir(intdir)
    end
  
    tot=min(nd,i+skips); %fewer pairs for later dates, no more than total # dates
    for j=i+1:tot
        corfile_small=[cordir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
        intfile_small=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int'];
         wgtfile=['wgtdir/' dates(i).name '_' dates(j).name '.wgt'];
        if(~exist(wgtfile,'file'))
            command=['imageMath.py -e=''arg(a*conj(b)*conj(c)*f)'' -o ' wgtfile ' --a=''' dates(i).slc ''' --b=''' dates(j).slc ''' --c=''' dates(i).slcVH ''' --f=''' dates(j).slcVH ''''];
            system(command);
        end
        
    end
end



for i=1:nd-1
      cordir=(['cordir' pol '/' dates(i).name '/']);
    intdir=(['intdir' pol '/' dates(i).name '/']);
  
    tot=min(nd,i+skips); %fewer pairs for later dates, no more than total # dates
    for j=i+1:tot
        corfile_small=[cordir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
        intfile_small=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int'];
         wgtfile=['wgtdir/' dates(i).name '_' dates(j).name '.wgt'];
        wgtfile='wgtdir/average.amp';
        if(~exist(corfile_small,'file'))
            disp(['running ' corfile_small])
            make_intcor_downlook(dates(i).slc,dates(j).slc,corfile_small,intfile_small,nx,ny,rlooks,alooks,1,wgtfile)
        else
            disp([corfile_small ' already made'])
        end
             
    end
end

if(1)
    %TEMPORARY mask based on downlooked geoms
    geomfile='merged/geom_master/hgt.rdr.1alks_3rlks.full';
    for i=1:nd-1
        cordir=(['cordir' pol '/' dates(i).name '/']);
        intdir=(['intdir' pol '/' dates(i).name '/']);
        tot=min(nd,i+skips); %fewer pairs for later dates, no more than total # dates
        for j=i+1:tot
            corfile_small=[cordir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
            intfile_small=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int']
            fidm=fopen(geomfile,'r');
            fid1=fopen(intfile_small,'r');
            fid2=fopen(corfile_small,'r');
            fido1=fopen('tmp1','w');
            fido2=fopen('tmp2','w');
            for k=1:newny
                a=fread(fidm,newnx,'real*8');
                b=fread(fid1,newnx,'real*4');
                c=fread(fid2,newnx,'real*4');
                b(a<0)=0;
                c(a<0)=0;
                fwrite(fido1,b,'real*4');
                fwrite(fido2,c,'real*4');
            end
            fclose('all');
            movefile('tmp1',intfile_small);
            movefile('tmp2',corfile_small);
        end
    end
end
%average cors
avgcorfile='cordir_VV/average.cor';
fido=fopen(avgcorfile,'w');
for j=1:newny
    fwrite(fido,zeros(newnx,1),'real*4');
end
fclose(fido)

for i=1:nd-1
    i
    j=i+1;
    cordir=(['cordir' pol '/' dates(i).name '/']);
    corfile_small=[cordir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
    fid1=fopen(avgcorfile,'r');
    fid2=fopen(corfile_small,'r');
    fido=fopen('tmpcor','w');
    for j=1:newny
        tmp1=fread(fid1,newnx,'real*4');
        tmp2=fread(fid2,newnx,'real*4');
        fwrite(fido,tmp1+tmp2/(nd-1),'real*4');
    end
    fclose('all');
    movefile('tmpcor',avgcorfile);
end

%average wgts
avgwgtfile='wgtdir/average.wgt';
fido=fopen(avgwgtfile,'w');
for j=1:ny
    fwrite(fido,zeros(nx*2,1),'real*4');
end
fclose(fido);

im=sqrt(-1);
for i=1:nd-1
    i
    j=i+1;
    wgtfile=['wgtdir/' dates(i).name '_' dates(j).name '.wgt'];

    fid1=fopen(avgwgtfile,'r');
    fid2=fopen(wgtfile,'r');
    fido=fopen('tmpwgt','w');
    for j=1:ny
        tmp1=fread(fid1,nx*2,'real*4');
        tmp2=fread(fid2,nx,'real*4');
        tmp1=tmp1(1:2:end)+im*tmp1(2:2:end);
        tmp2=exp(im*tmp2)/(nd-1);
        tmp2=tmp1+tmp2;
        tmp1=zeros(nx*2,1);
        tmp1(1:2:end)=real(tmp2);
        tmp1(2:2:end)=imag(tmp2);
        fwrite(fido,tmp1,'real*4');
    end
    fclose('all');
    movefile('tmpwgt',avgwgtfile);
end


return


for i=1:nd-1
    cordir=(['cordir' pol '/' dates(i).name '/']);
    intdir=(['intdir' pol '/' dates(i).name '/']);
    tot=min(nd,i+skips); %fewer pairs for later dates, no more than total # dates
    
    for j=i+1:tot
        intfile_small=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int']; 
        intfile_filt=[intdir dates(i).name '_' dates(j).name '_filt_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int'];
        intfile_2pi=[intdir dates(i).name '_' dates(j).name '_2pi_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
        intfile_filtunw=[intdir dates(i).name '_' dates(j).name '_filt_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
        intfile_unw=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
        if(~exist(intfile_unw,'file'))
            %filter
            command=['imageMath.py -e=''exp(I*a)'' -t cfloat -o tmp --a=''' intfile_small  ';' num2str(newnx) ';float;1;BSQ'''];
            system(command);
            command=['remove_nan tmp tmp2 ' num2str(newnx) ' ' num2str(newny) ' > /dev/null'];
            system(command);
            movefile('tmp2','tmp');
            command=['psfilt tmp ' intfile_filt ' ' num2str(newnx)];
            system(command);
            %unwrap filtered
            chdir('snaphu')
            command=['psfilt ../' intfile_filt ' snaphu.in ' num2str(newnx)]; %filtfilt
            system(command);
            command=['snaphu -f snaphu.conf'];
            system(command);
            command=['imageMath.py -e=''round((b-arg(a))/2/PI)'' -t short -n -o snaphu.2pi --a=''snaphu.in;' num2str(newnx) ';cfloat;1;BSQ'' --b=''snaphu.out;' num2str(newnx) ';float;1;BSQ'''];
            system(command);
            movefile('snaphu.2pi',['../' intfile_2pi]);
            chdir('..');
            %add 2pis to unfiltered
            command=['imageMath.py -e=''a+2*PI*b'' -o tmpunw -t float --a=''' intfile_small ';' num2str(newnx) ';float;1;BSQ'' --b=''' intfile_2pi  ';' num2str(newnx) ';short;1;BSQ'''];
            system(command);
            command=['imageMath.py -e=''a-round((a-b)/2/PI)*2*PI'' -o ' intfile_unw ' -t float --a=tmpunw --b=''snaphu/snaphu.out;' num2str(newnx) ';float;1;BSQ'''];
            system(command);
        end

    end
end