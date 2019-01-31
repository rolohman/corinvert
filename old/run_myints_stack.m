%parpool(8)
decide_ints_stack
slcdir='merged/SLC/';

for i=1:nd
    dates(i).name    = files(i).name(1:8);    
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
end

if(~exist('intdir'))
    mkdir('intdir')
end
rlooks=20;
alooks=8;
for i=1:nd-1
    for j=i+1
        intdir=(['intdir/' dates(i).name '/']);
        if(~exist(intdir,'dir'))
            mkdir(intdir)
        end
      
        intfile_small=[intdir dates(i).name '_' dates(j).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.int'];
        if(~exist(intfile_small,'file'))
            disp(['running ' intfile_small])         
            make_int_downlook(dates(i).slc,dates(j).slc,intfile_small,nx,ny,rlooks,alooks,0)
        else
            disp([intfile_small ' already made'])
        end
    end
end

