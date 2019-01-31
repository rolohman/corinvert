%parpool(8)
decide_ints_stack
slcdir='merged/SLC/';

for i=1:nd
    dates(i).name    = files(i).name(1:8);    
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc'];
end

if(~exist('cordir'))
    mkdir('cordir')
end

for i=1:nd-1
    for j=i+1:nd
        cordir=(['cordir/' dates(i).name '/']);
        if(~exist(cordir,'dir'))
            mkdir(cordir)
        end
        intcordir=[cordir dates(j).name '/'];
        if(~exist(intcordir,'dir'))
            mkdir(intcordir);
        end
        corfile_small=[intcordir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
        if(~exist(corfile_small,'file'))
            disp(['running ' corfile_small])         
            mycor_slcs_downlook(dates(i).slc,dates(j).slc,corfile_small,nx,ny,rlooks,alooks,0)
        else
            disp([corfile_small ' already made'])
        end
    end
end

