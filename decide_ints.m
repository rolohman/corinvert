home=pwd;
if(regexp(home,'156'))
    pth=156;
    nx=50206;%swaths 2&3
    ny=13123;
elseif(regexp(home,'54'))
    pth=54;
    nx=22502;%swath 1
    ny=12229;
elseif(regexp(home,'47'))
    pth=47;
    nx=47355;%swaths 1&2
    ny0=12702;
    ny=24780;
elseif(regexp(home,'166')) %mojave
    pth=166;
    nx=68912;
    ny=7282;
elseif(regexp(home,'143')) %Houston
    pth=143;
    nx=67522;
    ny=14511;
    masterdate='20170824';
elseif(regexp(home,'115')) %Yakima
    pth=115;
    nx=21219;
    ny=4179;
    masterdate='20170518';
elseif(regexp(home,'42')); %Yakima
    pth=42;
    nx=24429;
    ny=5542
    masterdate='20170531';
elseif(regexp(home,'64')); %Yakima
    pth=64;
    nx=26016;
    ny=5535;
    masterdate='20170521';
else
    disp('run from path directories')
end
files=dir('*_tops.xml');
for i=1:length(files)
    dates(i).name=files(i).name(1:8);
    dates(i).dn=datenum(dates(i).name,'yyyymmdd');
end
if(pth==47)
    disp('not using last date for path 47')
    dates=dates(1:end-1);
end


dn=[dates.dn];

nd=length(dn);
ints=[];
if(exist('masterdate','var'))
    dn_master = datenum(masterdate,'yyyymmdd');
    testid    = find(dn==dn_master);
    if(length(testid)==1)
        j     = testid;
    else
        disp(['no date found matching ' masterdate]);
    end
else
    j=1;
end
testid=j
disp(['using master date=' dates(j).name])
  
for i=[1:testid-1 testid+1:nd]
    ints(end+1).id1=testid;
    ints(end).id2=i;
end
id1=[ints.id1];
id2=[ints.id2];
ni=length(id1);


for i=1:nd
    dates(i).bp=0;
    %dates(i).bpstd=0;
    if(i~=testid) %master date
        filedir=[dates(testid).name '/int_' dates(testid).name '_' dates(i).name];
        file=[filedir '/isce.log'];
        if(exist(filedir,'dir'))
            [a,b]=grep('-s','Bperp',file);
            if(length(a)>0)
                tmp=regexp(b.match,'=\s*([^A-Za-Z]+)','tokens','once');
                for j=1:length(tmp)
                    jnk(j)=str2num(char(tmp{j}));
                end
                dates(i).bp=mean(jnk);
                %dates(i).bpstd=std(jnk);
            else
                disp(['bp not found in ' file]);
            end
        else
            disp([dates(testid).name '-' dates(i).name ' not processed yet']);
        end
    end   
end
bp=[dates.bp];

