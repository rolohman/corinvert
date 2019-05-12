%parpool(4)
pol='_VV';
decide_ints_stack
px=7;
py=3;
im=sqrt(-1);

for i=1:nd
    dates(i).name    = files(i).name(1:8);
end

clear ints
for i=1:nd-1
    j=i+1;
    intdir_orig = (['intdir' pol '/' dates(i).name '/']);
    cordir_orig = (['cordir' pol '/' dates(i).name '/']);
    intdir = (['intdir2' pol '/' dates(i).name '/']);

    if(~exist(intdir,'dir'))
        mkdir(intdir)
    end
    name   = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).orig  = [intdir_orig name '_highpass2_fix.unw']; %unfiltered unwrapped
    ints(i).cor         = [cordir_orig name '.cor'];

    ints(i).msk         = [intdir name '.msk'];
    ints(i).filt1   = ['smallfilt.int'];  %masked infilled with filtered
    ints(i).filtw1  = ['smallfiltw.int'];  %masked infilled with filtered
    ints(i).filt2   = ['bigfilt.int']; %psfilt version for unwrapping
    ints(i).psfilt  = [intdir name '_psfilt.int']; %psfilt version for unwrapping
    ints(i).pis     = [intdir name '_2pi.unw'];
    ints(i).unw     = [intdir name '.unw']; %unfiltered unwrapped
    ints(i).long    = [intdir name '_low.unw'];  %long-wavelength component
    ints(i).deramp  = [intdir name '_highpass.unw']; %unfiltered unwrapped
end



wgtfile=['intdir' pol '/average.' num2str(alooks) 'alks_' num2str(rlooks) 'rlks.cor'];
wgtfile=['longfiltmask.i1']
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
        mask1     = fread(fid1,[newnx newny],'integer*1');
        mask2     = fread(fid2,[newnx newny],'real*4');
        mask      = and(mask1,mask2>0.35);   
        
        fid1      = fopen(ints(i).msk,'w');
        fid2      = fopen('snaphu/snaphu.msk','w');
        fwrite(fid1,mask2,'integer*1');
        fwrite(fid2,mask,'real*4');
        fclose('all');
        clear mask mask1 mask2
        
        %filter twice and fill masked area
        myfilt(ints(i).orig,ints(i).msk,ints(i).filt1,3,3,newnx,newny,2,1,1,ints(i).filtw1);
        myfilt(ints(i).orig,ints(i).msk,ints(i).filt2,45,45,newnx,newny,2,1,1,'/dev/null');
        
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

%remove long-wavelength signals
for i=1:nd-1
    if(~exist(ints(i).long2,'file'))
        disp(['need to run ' ints(i).long2])
        %filter long wavelengths
        myfilt(ints(i).unw,'longfiltmask.i1',ints(i).long2,200,200,newnx,newny,2,3,4,'/dev/null');
        %myfilt(ints(i).unw,ints(i).msk,ints(i).long,200,200,newnx,newny,2,3,4,'/dev/null');
        
        %remove long wavelength from unw
        command=['imageMath.py -e=''a-b'' -t float -n -o ' ints(i).deramp2 ' --a=''' ints(i).unw  ';' num2str(newnx) ';float;1;BSQ'' --b=''' ints(i).long2  ';' num2str(newnx) ';float;1;BSQ'''];
        system(command);
    else
        disp(['already highpass filtered ' ints(i).int])
    end
end
