%parpool(8)
pol='_VV';
decide_ints_stack
slcdir=['merged/SLC' pol '/'];
skips=2;  %1=sequential, larger=longer time pairs
for i=1:nd
    dates(i).name    = files(i).name(1:8);    
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
end

if(~exist(['cordir' pol]))
    mkdir(['cordir' pol])
end
if(~exist(['intdir' pol]))
    mkdir(['intdir' pol])
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
    if(i==1)
        tot=nd;  %do all pairs vs. master date
    else
        tot=min(nd,i+skips); %fewer pairs for later dates, no more than total # dates
    end
    for j=i+1:tot
        corfile_small=[cordir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
        intfile_small=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int'];
        if(~exist(corfile_small,'file'))
            disp(['running ' corfile_small])
            make_intcor_downlook(dates(i).slc,dates(j).slc,corfile_small,intfile_small,nx,ny,rlooks,alooks,0)
        else
            disp([corfile_small ' already made'])
        end
    end
end

