function download_orbits(dn,dat_ab,sentime)
%first look for precise orbits
orbitdir='/data/Sentinel/precise/';

dat_sat(dat_ab==1)='A';
dat_sat(dat_ab==2)='B';


for iters=1:2 %do not remember why we iterate twice.
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
        tmp=find(and(and(odn1<dn(i),odn2>dn(i)),dat_ab(i)==sab));
        foundorbits(i)=length(tmp);
    end
    
    bado=find(foundorbits==0);
    tmp=urlread('https://s1qc.asf.alaska.edu/aux_poeorb/');
    for i=1:length(bado)
        a=regexp(tmp,['S1' dat_sat(bado(i)) '.{38}V' datestr(dn(bado(i))-1,'yyyymmdd')]);
        if(a)
            filename=tmp(a(1):a(1)+76)
            path=['https://s1qc.asf.alaska.edu/aux_poeorb/' filename];
            system(['wget ' path ' --no-check-certificate -O ' orbitdir filename]);
        else
            disp(['not found in poe:' datestr(dn(bado(i)))]);
        end
    end
end

bado=find(foundorbits==0);
tmp=urlread('https://s1qc.asf.alaska.edu/aux_resorb/');
for i=1:length(bado)
    test=['S1' dat_sat(bado(i)) '.{38}V' datestr(dn(bado(i)),'yyyymmdd')];
    a=regexp(tmp,test);
    if(a)
        
        for j=1:length(a)
            filename=tmp(a(j):a(j)+76);
            t1=str2num(filename(52:57));
            t2=str2num(filename(68:73));
            if(and(sentime(bado(i))>t1,sentime(bado(i))<t2))
                disp(['found! ' filename])
                path=['https://s1qc.asf.alaska.edu/aux_resorb/' filename];
                system(['wget ' path ' --no-check-certificate -O ' orbitdir filename]);
                break
            end
        end
        
    else
        disp(['not found in res:' datestr(dn(bado(i)))]);
    end
end