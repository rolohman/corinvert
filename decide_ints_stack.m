home=pwd;
params

if(~exist('pol','var'))
    pol='';
end
slcdir=['merged/SLC' pol '/'];
files=dir([slcdir '2*']);
dates=[];
for i=1:length(files)
    dates(i).name=files(i).name(1:8);
    dates(i).dn=datenum(dates(i).name,'yyyymmdd');
end

nd     = length(dates);
ints=[];
if(nd==0)
    
else
    for j=1:nd-1
        for i=j+1:nd %make all
            ints(end+1).id1=j;
            ints(end).id2=i;
        end
    end
    id1    = [ints.id1]';
    id2    = [ints.id2]';
    diags  = find(id2==id1+1);
    ni     = length(ints);
    
    dn     = [dates.dn];
    dn1    = dn(id1);
    dn2    = dn(id2);
    intdt  = diff(dn); %intervals between dates (not interferograms)
end
rlooks = 15;
alooks = 4;
%rlooks=7;
%alooks=3;
%rlooks=3;
%alooks=1;
n      = rlooks*alooks;



if(exist('nx','var'))
    newnx  = floor(nx/rlooks)
    newny  = floor(ny/alooks);
else
    file=['merged/SLC' pol '/' dates(1).name '/' dates(1).name '.slc.full.vrt'];
    if(exist(file))
        [a,b]=system(['grep rasterXSize ' file]);
        tmp=regexp(b,'rasterXSize="(\d+)" rasterYSize="(\d+)">','tokens');
        if(length(tmp)==1)
            nx=str2num(tmp{1}{1});
            ny=str2num(tmp{1}{2});
            newnx  = floor(nx/rlooks);
            newny  = floor(ny/alooks);
            system(['echo ''nx=' num2str(nx) ';'' >>params.m']);
            system(['echo ''ny=' num2str(ny) ';'' >>params.m']);
            system(['echo ''newnx=' num2str(newnx) ';'' >>params.m']);
            system(['echo ''newny=' num2str(newny) ';'' >>params.m']);
            
        else
            disp(['no match for width/length in: ' b]);
        end
    else
        disp('no width/length info found yet, first interferogram not made');
    end
end
