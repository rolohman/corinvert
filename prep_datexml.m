%path='/data/rlohman/Sentinel/Chile/47_1096/data/';
%path='/data/rlohman/Sentinel/Chile/156_674/data/';
%path='/data/rlohman/Sentinel/Chile/54_673/data/';
%swath=2;
home=pwd;
path=[home '/data/'];
npath=length(path);
system(['ls -d1 ' path 'S*zip > SAFE.txt'])
%system('sed -i -e ''s/zip/SAFE/g'' SAFE.txt')
cskh5=textread('SAFE.txt','%s');
cskh5=cellstr(cskh5);

%sort SAFE file names, otherwise SDV -SSV modes will mess up the files

for i=1:length(cskh5)   
    hh=cskh5{i};
   hh=hh(npath+1:end); %remove root path
   sentime(i)=str2num(hh(27:32));

%   sentime(i)=str2num(hh(69:74));
    if(regexp(hh,'S1A'))
        datsab(i)=1;
        datsbs(i)='A';
    elseif(regexp(hh,'S1B'))
        datsab(i)=2;
        datsbs(i)='B';
    end
    hh=regexp(cskh5{i},'SLC__...._(\d{8})','tokens');
    if(length(hh)==1)
        dnall(i)=datenum(char(hh{1}),'yyyymmdd');
    else
        disp('bad number of date matches')
        return
    end
  end
[dn,uid]=unique(dnall);
datsab=datsab(uid);
datsbs=datsbs(uid);
sentime=sentime(uid);

for iters=1:2
orbits=dir('/data/Sentinel/precise/S*');
numo=length(orbits);
for i=1:numo
    name=orbits(i).name;
    odn1(i)=datenum(name(43:50),'yyyymmdd');
    odn2(i)=datenum(name(59:66),'yyyymmdd');
    tmp=name(3);
    switch tmp
        case 'A'
            sab(i)=1;
        case 'B'
            sab(i)=2;
    end
end



foundorbits=zeros(1,length(dn));
for i=1:length(dn)
    tmp=find(and(and(odn1<dn(i),odn2>dn(i)),datsab(i)==sab));
foundorbits(i)=length(tmp);
end
    
bado=find(foundorbits==0);
tmp=urlread('https://s1qc.asf.alaska.edu/aux_poeorb/');
for i=1:length(bado)
    a=regexp(tmp,['S1' datsbs(bado(i)) '.{38}V' datestr(dn(bado(i))-1,'yyyymmdd')]);
    if(a)
        filename=tmp(a(1):a(1)+76)
        path=['https://s1qc.asf.alaska.edu/aux_poeorb/' filename];
        system(['wget ' path ' --no-check-certificate -O /data/Sentinel/precise/' filename]);
    else
        disp(['not found in poe:' datestr(dn(bado(i)))]);
    end
end
end

bado=find(foundorbits==0);
tmp=urlread('https://s1qc.asf.alaska.edu/aux_resorb/');
for i=1:length(bado)
    test=['S1' datsbs(bado(i)) '.{38}V' datestr(dn(bado(i)),'yyyymmdd')];
    a=regexp(tmp,test);
    if(a)
        
        for j=1:length(a)
            filename=tmp(a(j):a(j)+76);
            t1=str2num(filename(52:57));
            t2=str2num(filename(68:73));
            if(and(sentime(bado(i))>t1,sentime(bado(i))<t2))
                disp(['found! ' filename])
                path=['https://s1qc.asf.alaska.edu/aux_resorb/' filename];
                system(['wget ' path ' --no-check-certificate -O /data/Sentinel/precise/' filename]);
                break
            end
        end
        
    else
        disp(['not found in poe:' datestr(dn(bado(i)))]);
    end
end

for i=1:length(dn)
    date=datestr(dn(i),'yyyymmdd');
    fid3 = fopen([date,'_tops.xml'],'w');
    id=find(dnall==dn(i));
    
    if length(id)==1;
        files=cskh5{id};
    else
        files=['"' strjoin(cskh5(id),'","') '"'];
    end
    fprintf(fid3,'<?xml version="1.0" encoding="UTF-8"?> \n');
    fprintf(fid3,'<component name="master"> \n');
    fprintf(fid3,'  <property name="safe">%s</property>\n',files);
%    fprintf(fid3,'  <property name="swath number">%d</property>\n',swath);
    fprintf(fid3,'  <property name="output directory">%s</property>\n',date);
    fprintf(fid3,'  <property name="orbit directory">/data/Sentinel/precise</property>\n');
    fprintf(fid3,'  <property name="auxiliary data directory">/data/Sentinel/aux_cal</property>\n');
    fprintf(fid3,'</component>\n');
    fclose(fid3);
end


