%parpool(4)
pol='_VV';
decide_ints_stack
px=5;
py=3;
im=sqrt(-1);
slcdir=['merged/SLC' pol '/'];


if(~exist(['wgtdir']))
    mkdir(['wgtdir'])
end
if(~exist(['intdir' pol]))
    mkdir(['intdir' pol])
end


for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
end

slcdirVH=['merged/SLC_VH/'];
for i=1:nd
    dates(i).slcVH     = [slcdirVH dates(i).name '/' dates(i).name '.slc.full'];
end


clear ints
for i=1:nd-1
    j=i+1;
     intdir = (['intdir' pol '/' dates(i).name '/']);
     if(~exist(intdir,'dir'))
        mkdir(intdir)
    end
      name   = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).int         = [intdir name '.int'];
    ints(i).cor         = [intdir name '.cor'];
    ints(i).int2        = [intdir name '_' num2str(px) '.int'];
    ints(i).cor2        = [intdir name '_' num2str(px) '.cor'];
    ints(i).msk         = [intdir name '.msk'];
    ints(i).filt1   = ['smallfilt.int'];  %masked infilled with filtered
    ints(i).filtw1  = ['smallfiltw.int'];  %masked infilled with filtered
    ints(i).filt2   = ['bigfilt.int']; %psfilt version for unwrapping
    ints(i).psfilt  = [intdir name '_psfilt.int']; %psfilt version for unwrapping
    ints(i).pis     = [intdir name '_2pi.unw'];
    ints(i).unw     = [intdir name '.unw']; %unfiltered unwrapped
    ints(i).long    = [intdir name '_low.unw'];  %long-wavelength component
    ints(i).deramp  = [intdir name '_highpass.unw']; %unfiltered unwrapped
    ints(i).vvvh    = ['wgtdir/' dates(i).name '_' dates(j).name '.int'];

end

wgtfile=['intdir' pol '/average.cor'];
geomfile='merged/geom_master/hgt.rdr.1alks_3rlks.full';

%VVVH stuff
if(~exist(wgtfile,'file'))
    for i=1:nd-1
        j=i+1;
        command=['imageMath.py -e=''arg(a*conj(b)*conj(c)*f)'' -o ' ints(i).vvvh ' --a=''' dates(i).slc ''' --b=''' dates(j).slc ''' --c=''' dates(i).slcVH ''' --f=''' dates(j).slcVH ''''];
        system(command);
    end
    
    fido=fopen(wgtfile,'w');
    for i=1:nd-1
        fid(i)=fopen(ints(i).vvvh,'r');
    end
    vvvh=nan(nd-1,nx);
    for j=1:ny
        for i=1:nd-1
            vvvh(i,:) = fread(fid(i),nx,'real*4');    
        end
        vvvh=exp(im*vvvh);
        fwrite(fido,abs(mean(vvvh,1,'omitnan')),'real*4');
    end
    fclose('all');    
end
copyxml(ints(i).vvvh,wgtfile);
command=['looks.py -i ' wgtfile ' -r ' num2str(rlooks) ' -a ' num2str(alooks)];
system(command);



wgtfile=['intdir' pol '/average.cor'];
geomfile='merged/geom_master/hgt.rdr.1alks_3rlks.full';

%make ints and corfiles if not already made
for i=1:nd-1
    j=i+1;
    if(~exist(ints(i).int,'file'))
        disp(['running ' ints(i).int])
        make_intcor_downlook(dates(i).slc,dates(j).slc,ints(i).cor,ints(i).int,nx,ny,rlooks,alooks,1,wgtfile)
        make_intcor_anyposting(dates(i).slc,dates(j).slc,ints(i).cor2,ints(i).int2,nx,ny,rlooks,alooks,px,py,1,wgtfile)
        
        mask_ztopo(geomfile,ints(i).int,1,newnx,newny)
        mask_ztopo(geomfile,ints(i).int2,1,newnx,newny)
        mask_ztopo(geomfile,ints(i).cor,1,newnx,newny)
        mask_ztopo(geomfile,ints(i).cor2,1,newnx,newny)

    else
        disp([ints(i).int ' already made'])
    end
end

wgtfile=['intdir' pol '/average.' num2str(alooks) 'alks_' num2str(rlooks) 'rlks.cor'];
%now the filtering/unwrapping section
for i=1:nd-1
    j=i+1;
    
    if(~exist(ints(i).unw,'file'))
        disp(['unwrapping ' ints(i).int])
        delete('maskfill.int'); %make sure this is blank
        delete('snaphu/snaphu.out');
        delete('snaphu/snaphu.in');
        delete('snaphu/snaphu.msk');
        %make mask
        fid1      = fopen(wgtfile,'r');
        fid2      = fopen(ints(i).cor2,'r');
        mask1     = fread(fid1,[newnx newny],'real*4');
        mask2     = fread(fid2,[newnx newny],'real*4');
        mask      = and(mask1>0.1,mask2>0.35);   
        
        fid1      = fopen(ints(i).msk,'w');
        fid2      = fopen('snaphu/snaphu.msk','w');
        fwrite(fid1,mask2,'integer*1');
        fwrite(fid2,mask,'real*4');
        fclose('all');
        clear mask mask1 mask2
        
        %filter twice and fill masked area
        myfilt(ints(i).int,ints(i).msk,ints(i).filt1,3,3,newnx,newny,2,1,1,ints(i).filtw1);
        myfilt(ints(i).int,ints(i).msk,ints(i).filt2,40,40,newnx,newny,2,1,1,'/dev/null');
        
        mask_ztopo(geomfile,ints(i).filt1,1,newnx,newny)
        mask_ztopo(geomfile,ints(i).filt2,1,newnx,newny)

        fid0=fopen(ints(i).msk,'r');
        fid1=fopen(ints(i).int,'r');
        fid2=fopen(ints(i).filt1,'r');
        fid3=fopen(ints(i).filtw1,'r');
        fid4=fopen(ints(i).filt2,'r');
        fido=fopen('maskfill.int','w');
        
        outc          = zeros(1,newnx*2);
        for k=1:newny
            m  = fread(fid0,newnx,'integer*1');
            a  = fread(fid1,newnx,'real*4');
            b  = fread(fid2,newnx,'real*4');
            c  = fread(fid3,newnx,'real*4');
            d  = fread(fid4,newnx,'real*4');
            
            b  = exp(im*b);
            d  = exp(im*d);
            m1 = m==1;
            
            %a, unless masked, then b weighted by distance+d, unless
            %masked, then d.
            out(m1)           = a(m1);
            out(~m1)          = angle(b(~m1).*c(~m1)+d(~m1).*(1-c(~m1)));
            out(isnan(b))     = angle(d(isnan(b)));
            out(isnan(a))     = NaN;
            outc(1:2:end)     = cos(out);
            outc(2:2:end)     = sin(out);
            outc(isnan(outc)) = 0;
            fwrite(fido,outc,'real*4');
        end
        
        
        %psfilt, remove nans
        command=['psfilt maskfill.int ' ints(i).psfilt ' ' num2str(newnx) ' 0.25'];
        system(command);
        
        %unwrap filtered with snaphu
        copyfile(ints(i).psfilt,'snaphu/snaphu.in')
        %convert back to r4 to save space
        command=['imageMath.py -e=''arg(a)'' -t float -n -o ' ints(i).psfilt ' --a=''snaphu/snaphu.in;' num2str(newnx) ';cfloat;1;BSQ'''];
        system(command);
        
        chdir('snaphu')
        
        command=['snaphu -f snaphu.conf'];
        system(command);
        chdir('..')
        
        %add 2pis to unfiltered
        
        command=['imageMath.py -e=''round((b-a)/2/PI)'' -t short -n -o ' ints(i).pis ' --a=''' ints(i).int ';' num2str(newnx) ';float;1;BSQ'' --b=''snaphu/snaphu.out;' num2str(newnx) ';float;1;BSQ'''];
        system(command);
        command=['imageMath.py -e=''a+2*PI*b'' -o ' ints(i).unw ' -n -t float --a=''' ints(i).int ';' num2str(newnx) ';float;1;BSQ'' --b=''' ints(i).pis  ';' num2str(newnx) ';short;1;BSQ'''];
        system(command);
    else
        disp(['done unwrapping ' ints(i).int])
    end
end

%temporary addition, fixing long wavelength mask
% fid1=fopen('demerr.r4','r');
% fid2=fopen('sigstd.r4','r');
% fido=fopen('longfiltmask.i1','w');
% a=fread(fid1,[newnx newny],'real*4');
% b=fread(fid2,[newnx newny],'real*4');
% c=and(abs(a)<0.002,b<1);
% fwrite(fido,c,'integer*1');
% fclose('all');

%remove long-wavelength signals
for i=1:nd-1
    if(~exist(ints(i).long,'file'))
        disp(['need to run ' ints(i).long])
        %filter long wavelengths
        %myfilt(ints(i).unw,'longfiltmask.i1',ints(i).long,200,200,newnx,newny,2,3,4,'/dev/null');
        myfilt(ints(i).unw,ints(i).msk,ints(i).long,200,200,newnx,newny,2,3,4,'/dev/null');
        
        %remove long wavelength from unw
        command=['imageMath.py -e=''a-b'' -t float -n -o ' ints(i).deramp ' --a=''' ints(i).unw  ';' num2str(newnx) ';float;1;BSQ'' --b=''' ints(i).long  ';' num2str(newnx) ';float;1;BSQ'''];
        system(command);
    else
        disp(['already highpass filtered ' ints(i).int])
    end
end

% if(1)
%     %TEMPORARY mask based on downlooked geoms
%     geomfile='merged/geom_master/hgt.rdr.1alks_3rlks.full';
%     for i=1:nd-1
%         cordir=(['cordir' pol '/' dates(i).name '/']);
%         intdir=(['intdir' pol '/' dates(i).name '/']);
%         tot=min(nd,i+skips); %fewer pairs for later dates, no more than total # dates
%         for j=i+1:tot
%             intfile_unw=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
%             fidm=fopen(geomfile,'r');
%             fid1=fopen(intfile_unw,'r');
%             fido1=fopen('tmp1','w');
%             for k=1:newny
%                 a=fread(fidm,newnx,'real*8');
%                 b=fread(fid1,newnx,'real*4');
%                 b(a<0)=0;
%                 fwrite(fido1,b,'real*4');
%             end
%             fclose('all');
%             movefile('tmp1',intfile_unw);
%         end
%     end
% end
%
%
