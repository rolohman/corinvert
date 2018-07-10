home=pwd;
params

files=dir('datexml/*_tops.xml');
for i=1:length(files)
    dates(i).name=files(i).name(1:8);
    dates(i).dn=datenum(dates(i).name,'yyyymmdd');
end

nd     = length(dates);

ints=[];
for j=1:nd-1
    for i=j+1:nd
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
rlooks = 7;
alooks = 3;
n      = rlooks*alooks;

    

if(exist('nx','var'))
    newnx  = floor(nx/rlooks)
    newny  = floor(ny/alooks);
else
    file='slcs/master.slc.full.vrt';
    if(exist(file))
        [a,b]=system(['grep rasterXSize ' file]);
        tmp=regexp(b,'rasterXSize="(\d+)" rasterYSize="(\d+)">','tokens');
        if(length(tmp)==1)
            nx=str2num(tmp{1}{1});
            ny=str2num(tmp{1}{2});
            newnx  = floor(nx/rlooks)
            newny  = floor(ny/alooks);
            system(['echo ''nx=' num2str(nx) ';'' >>params.m']);
            system(['echo ''ny=' num2str(ny) ';'' >>params.m']);
            
        else
            disp(['no match for width/length in: ' b]);
        end
    else
        disp('no width/length info found yet, first interferogram not made');
    end
end