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
nd     = length(dates);
%assume all ints made
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
newnx  = floor(nx/rlooks)
newny  = floor(ny/alooks);

