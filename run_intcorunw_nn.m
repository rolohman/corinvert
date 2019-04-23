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
            command=['imageMath.py -e=''abs(arg(a*conj(b)*conj(c)*f))'' -o ' wgtfile ' --a=''' dates(i).slc ''' --b=''' dates(j).slc ''' --c=''' dates(i).slcVH ''' --f=''' dates(j).slcVH ''''];
            system(command);
        end
        if(~exist(corfile_small,'file'))
            disp(['running ' corfile_small])
            make_intcor_downlook(dates(i).slc,dates(j).slc,corfile_small,intfile_small,nx,ny,rlooks,alooks,1,wgtfile)
        else
            disp([corfile_small ' already made'])
        end
             
    end
end
return
for i=1:nd-1
    cordir=(['cordir' pol '/' dates(i).name '/']);
    intdir=(['intdir' pol '/' dates(i).name '/']);
      tot=min(nd,i+skips); %fewer pairs for later dates, no more than total # dates
  
   for j=i+1:tot
        corfile_small=[cordir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
        intfile_small=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int'];
 
   intfile_filt=[intdir dates(i).name '_' dates(j).name '_filt_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int'];
        intfile_unw=[intdir dates(i).name '_' dates(j).name '_filt_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.unw'];
       %filter
        command=['imageMath.py -e=''exp(j*a)'' -t cfloat -o tmp --a=''' intfile_small ''''];
        system(command);
        command=['psfilt tmp ' intfile_filt ' ' num2str(newnx)];
        system(command);
        %unwrap filtered
        command(['snaphu ' intfile_filt ' ' num2str(newnx) ' -s -o tmp.unw --tile 10 10 100 100']);
        system(command)
   end
end