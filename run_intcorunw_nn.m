%parpool(8)
pol='_VV';
decide_ints_stack
px=7;
py=3;
im=sqrt(-1);
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
    cordir2=(['cordir2' pol '/' dates(i).name '/']);
    intdir2=(['intdir2' pol '/' dates(i).name '/']);
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
        corfile2_small=[cordir2 dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
        intfile2_small=[intdir2 dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int'];
        %wgtfile=['wgtdir/' dates(i).name '_' dates(j).name '.wgt'];
        wgtfile='wgtdir/average.amp';
        if(~exist(corfile_small,'file'))
            disp(['running ' corfile_small])
            make_intcor_downlook(dates(i).slc,dates(j).sl,corfile_small,intfile_small,nx,ny,rlooks,alooks,1,wgtfile)
            make_intcor_anyposting(dates(i).slc,dates(j).slc,corfile2_small,intfile2_mapp,nx,ny,px,py,rlooks,alooks,1,wgtfile)           
            for x=1:1
                geomfile='merged/geom_master/hgt.rdr.1alks_3rlks.full';
                fidm=fopen(geomfile,'r');
                fid1=fopen(intfile_small,'r');
                fid2=fopen(corfile_small,'r');
                fid1=fopen(intfile2_small,'r');
                fid2=fopen(corfile2_small,'r');
                fido1=fopen('tmp1','w');
                fido2=fopen('tmp2','w');
                fido3=fopen('tmp3','w');
                fido4=fopen('tmp4','w');
                for k=1:newny
                    a=fread(fidm,newnx,'real*8');
                    b=fread(fid1,newnx,'real*4');
                    c=fread(fid2,newnx,'real*4');
                    d=fread(fid1,newnx,'real*4');
                    e=fread(fid2,newnx,'real*4');
                    b(a<0)=0;
                    c(a<0)=0;
                    d(a<0)=0;
                    e(a<0)=0;
                    fwrite(fido1,b,'real*4');
                    fwrite(fido2,c,'real*4');
                    fwrite(fido3,d,'real*4');
                    fwrite(fido4,e,'real*4');
                end
                fclose('all');
                movefile('tmp1',intfile_small);
                movefile('tmp2',corfile_small);
                movefile('tmp3',intfile2_small);
                movefile('tmp4',corfile2_small);
                else
                    disp([corfile_small ' already made'])
            end
        end
    end
end

%average cors
avgcorfile='cordir_VV/average.cor';
if(~exist(avgcorfile,'file'))
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
end

%average wgts
avgwgtfile='wgtdir/average.wgt';
if(~exist(avgwgtfile,'file'))
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
end


for i=1:nd-1
    j=i+1;
    cordir = (['cordir2' pol '/' dates(i).name '/']);
    intdir = (['intdir' pol '/' dates(i).name '/']);
    name   = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    intfile         = [intdir name '.int'];
    corfile         = [cordir name '.cor'];
    intmask         = [intdir name '.msk'];
    intfile_long    = [intdir name 'lonw.int'];  %long-wavelength component
    intfile_filt1   = [intdir name 'filt.int'];  %masked infilled with filtered
    intfile_filt2   = [intdir name 'psfilt.int']; %psfilt version for unwrapping
    intfile_filtunw = [intdir name 'psfilt.unw']; %psfilt version for unwrapping

    intfile_2pi     = [intdir name '2pi.unw'];
    intfile_unw     = [intdir name '.unw']; %unfiltered unwrapped
    
    if(~exist(intfile_unw,'file'))
        %make mask
        fid       = fopen(avgwgtfile,'r');
        mask1     = fread(fid,[newnx newny],'real*4');
        fid       = fopen(corfile,'r');
        mask2     = fread(fid,[newnx newny],'real*4');
        mask      = and(mask1>0.15,mask2>0.4);
        fid       = fopen(intmask,'w');
        fwrite(fid,mask,'bit1');
        fclose('all');
        
        
        %filter long wavelengths
        

        ry=500;
        windx=ones(1,ry*2+1);
        windy=ones(1,ry*2+1);
        windx=windx/sum(windx);
        windy=windy/sum(windy);
        for i=1:nd-1
            i
            if(~exist(ints(i).filt,'file'))
                
                z       = zeros(1,newnx);
                fidi    = fopen(ints(i).unw,'r');
                fido    = fopen(ints(i).filt,'w');
                
                mskunw=fread(fidi,[newnx,ry],'real*4');
                mskunw(~mask(:,1:ry))=0; %mask on fly
                mskunw = [flipud(mskunw');nan(ry+1,newnx)];
                mskunw(isnan(mskunw))=0;
                
                %now go to end (passing end by ry lines, filling in with zeros. "active"
                %line is at ry+1th row
                
                for j=1:newny
                    
                    mskunw=circshift(mskunw,1);
                    [a,count1]=fread(fidi,newnx,'real*4');
                    a(isnan(a))=0;
                    if(count1==newnx)
                        a(~mask(:,ry+j))=0;
                        mskunw(1,:)=a;
                    else
                        mskunw(1,:)=z;
                    end
                    good = mskunw~=0;
                    a    = good;
                    c    = mskunw;
                    
                    a    = windy*a;
                    c    = windy*c;
                    
                    asum = conv(a,windx,'same');
                    csum = conv(c,windx,'same');
                    out  = csum./asum;
                    
                    fwrite(fido,out,'real*4'); %1000pixel filtered product, at all pixels, even masked ones.
                end
                fclose(fidi);
                fclose(fido);
                fidi1   = fopen(ints(i).unw,'r');
                fidi2   = fopen(ints(i).filt,'r');
                fido    = fopen(ints(i).fixunw,'w');
                for j=1:newny
                    a=fread(fidi1,newnx,'real*4');
                    b=fread(fidi2,newnx,'real*4');
                    fwrite(fido,a-b,'real*4');
                end
            end
        end

        
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








if(1)
    %TEMPORARY mask based on downlooked geoms
    geomfile='merged/geom_master/hgt.rdr.1alks_3rlks.full';
    for i=1:nd-1
        cordir=(['cordir' pol '/' dates(i).name '/']);
        intdir=(['intdir' pol '/' dates(i).name '/']);
        tot=min(nd,i+skips); %fewer pairs for later dates, no more than total # dates
        for j=i+1:tot
            intfile_unw=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
            fidm=fopen(geomfile,'r');
            fid1=fopen(intfile_unw,'r');
            fido1=fopen('tmp1','w');
            for k=1:newny
                a=fread(fidm,newnx,'real*8');
                b=fread(fid1,newnx,'real*4');
                b(a<0)=0;
                fwrite(fido1,b,'real*4');
            end
            fclose('all');
            movefile('tmp1',intfile_unw);
        end
    end
end




for i=1:nd-1
    j=i+1
    intfile_unw1=['tsdir/fixunw/' dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
    intfile_unw2=['tsdir/resnaphu/' dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
    if(~exist(intfile_unw2,'file'))
        %filter
        command=['remove_nan ' intfile_unw1 ' tsdir/resnaphu/snaphu.in ' num2str(newnx) ' ' num2str(newny) ' > /dev/null'];
        system(command);
        chdir('tsdir/resnaphu')
        command=['snaphu -f snaphu.conf'];
        system(command);
        movefile('snaphu.out',intfile_unw);
        chdir('../..');
        
    end
end

