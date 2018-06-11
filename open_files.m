if(~exist(adir,'dir'))
    disp(['creating directory ' adir]);
    mkdir(adir)
end

for i=1:ni
    cordir=(['cordir/' dates(id1(i)).name '/']);
    intcordir=[cordir dates(id2(i)).name '/'];
    corfile_small=[intcordir dates(id1(i)).name '_' dates(id2(i)).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
    if(exist(corfile_small,'file'))
        fidi(i)=fopen(corfile_small,'r');
    else
        disp(['where is ' corfile_small '?'])
    end
end

%find last line written
if(exist([adir 'c0.cor'],'file'))
    fid0=fopen([adir 'c0.cor'],'r');
    for j=1:newny
        [tmp,count]=fread(fid0,newnx,'real*4');
        if(count<newnx)
            break
        end
    end
    fclose(fid0);
    
    online=j-1
    if(online<1)
        disp('c0.cor empty?')
        return
    end
    fid0=fopen([adir 'c0.cor'],'a');
    fid1=fopen([adir 'rms.cor'],'a');
    fidmin=fopen([adir 'cmin.cor'],'w');
    
    for i=1:nd
        fidr(i)=fopen([adir 'rel_' dates(i).name '.cor'],'a');
    end
    for i=1:nd-1
        fidp(i)=fopen([adir 'perm_' dates(i).name '_' dates(i+1).name '.cor'],'a');
    end
    for i=1:ni
        fseek(fidi(i),online*newnx*4,-1);
    end
else
    online=0;
    for i=1:nd
        fidr(i)=fopen([adir 'rel_' dates(i).name '.cor'],'w');
    end
    for i=1:nd-1
        fidp(i)=fopen([adir 'perm_' dates(i).name '_' dates(i+1).name '.cor'],'w');
    end
    fid0=fopen([adir 'c0.cor'],'w');
    fid1=fopen([adir 'rms.cor'],'w');
    fidmin=fopen([adir 'cmin.cor'],'w');
end
