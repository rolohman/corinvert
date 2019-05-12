%parpool(4)
pol='_VV';
decide_ints_stack
px=7;
py=3;
im=sqrt(-1);


if(~exist(['intdir' pol '_iter']))
    mkdir(['intdir' pol '_iter'])
end

for i=1:nd
    dates(i).name    = files(i).name(1:8);
end

clear ints
for i=1:nd-1
    j=i+1;
    intdir_orig = (['intdir' pol '/' dates(i).name '/']);
    cordir_orig = (['cordir' pol '/' dates(i).name '/']);
    intdir = (['intdir' pol '/' dates(i).name '/']);

    if(~exist(intdir,'dir'))
        mkdir(intdir)
    end
    name   = [dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk'];
    ints(i).orig  = [intdir_orig name '_highpass_fix.unw']; %unfiltered unwrapped
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
     ints(i).deramp2  = [intdir name '_highpass2.unw']; %unfiltered unwrapped
    ints(i).long2    = [intdir name '_low2.unw'];  %long-wavelength component
 
end

%temporary addition, fixing long wavelength mask
fid1=fopen('dates_VV/simpledemerr.r4','r');
fid2=fopen('dates_VV/simpleres.r4','r');
fido=fopen('longfiltmask.i1','w');
a=fread(fid1,[newnx newny],'real*4');
b=fread(fid2,[newnx newny],'real*4');
c=and(abs(a)<0.001,b<1);
fwrite(fido,c,'integer*1');
fclose('all');

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
%
