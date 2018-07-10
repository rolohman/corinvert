home=pwd;
datapath=[home '/data/'];

datafiles=dir([datapath 'S*zip']);
nfiles=length(datafiles);

for i=1:nfiles
    name=datafiles(i).name;   
    if(regexp(name,'S1A'))
        dat_ab(i)=1;
    elseif(regexp(name,'S1B'))
        dat_ab(i)=2;
    end
    tmp=regexp(name,'SLC__...._(\d{8})T(\d{6})','tokens');
    if(length(tmp)==1)
        dnall(i)=datenum(tmp{1}{1},'yyyymmdd');
        sentime(i)=str2num(tmp{1}{2}); %utc
    else
        disp('bad number of date matches')
        return
    end
end
[dn,uid] = unique(dnall);
dat_ab   = dat_ab(uid);
sentime  = sentime(uid);


for i=1:length(dn)
    date=datestr(dn(i),'yyyymmdd');
    fid = fopen(['datexml/' date,'_tops.xml'],'w');
    id=find(dnall==dn(i));
    
    files='"';
    for j=1:length(id)-1
        files=[files datapath datafiles(j).name '","'];
    end
    files=[files datapath datafiles(j+1).name '"'];
    
    fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?> \n');
    fprintf(fid,'<component name="master"> \n');
    fprintf(fid,'  <property name="safe">%s</property>\n',files);
    fprintf(fid,'  <property name="output directory">%s</property>\n',date);
    fprintf(fid,'  <property name="orbit directory">/data/Sentinel/precise</property>\n');
    fprintf(fid,'  <property name="auxiliary data directory">/data/Sentinel/aux_cal</property>\n');
    fprintf(fid,'</component>\n');
    fclose(fid);
end

download_orbits(dn,dat_ab,sentime);
if(~exist('params.m','file'))
     system(['echo ''masterdate=''' datestr(dn(i),'yyyymmdd') ''';'' >params.m']);
end
    