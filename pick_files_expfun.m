function [dates,perms,c0s,dnr,fidi,fido,nx,ny]=pick_files_expfun(relDir,rdir,rlooks,alooks,pol)
%rflag = 1, open for writing/appending, 2, only for reading
dates = [];
perms = [];
c0s   = [];

if(~(rdir(end)=='/'))
    rdir=[rdir '/'];
end
if(or(rlooks>1,alooks>1))
    suff=['_' num2str(rlooks) 'r_' num2str(alooks) 'a'];
else
    suff='';
end
if(iscell(relDir))
	relDir=relDir{1};
end
tmp=dir([relDir '/rel*' suff '.cor.geo']);
for j=1:length(tmp)
    t=tmp(j).name;
    dates(end+1).name   = t(5:12);
    dates(end).filename = [tmp(j).folder '/' t];
    dates(end).dn       = datenum(t(5:12),'yyyymmdd');
    dates(end).origdir  = [tmp(j).folder];
end

tmp=dir([relDir '/perm*' suff '.cor.geo']);
for j=1:length(tmp)
    t=tmp(j).name;
    perms(end+1).name       = t(6:22);
    perms(end).filename     = [tmp(j).folder '/' t];
    perms(end).d1           = datenum(t(6:13),'yyyymmdd');
    perms(end).d2           = datenum(t(15:22),'yyyymmdd');
    perms(end).d            = mean([perms(end).d1 perms(end).d2]);
end

tmp=dir([relDir '/*c0' suff '.cor.geo']);
for j=1:length(tmp)
    t=tmp(j).name;
    c0s(end+1).filename     = [tmp(j).folder '/' t];
end
dn    = [dates.dn];
nd    = length(dn);

%get nx/ny from last file opened
vrt=[dates(end).filename '.vrt'];
if(exist(vrt,'file'))
    [a,b] = system(['grep rasterXSize ' vrt]);
   tmp   = regexp(b,'rasterXSize="(\d+)" rasterYSize="(\d+)"','tokens');
     if(length(tmp)==1)
        nx=str2num(tmp{1}{1}); ny=str2num(tmp{1}{2});
        disp(['width, length: ' num2str([nx ny])])
    end
else
    disp('no nx/ny info, need vrt file'),return
end

[jnk,sortid]=sort([dates.dn]);
dates = dates(sortid);

%load rain dates
if(exist('dates.list','file'))
    rdates=cellstr(num2str(load('dates.list')));
else
    disp('no date.list file, need dates of suspected rain events, YYYYMMDD')
    return
end

dnr    = datenum(rdates,'yyyymmdd');


goodr  = and(dnr<max(dn),dnr>min(dn)); %only use rain during time series
rdates  = rdates(goodr);
dnr    = dnr(goodr);
nr     = length(dnr);
  
%inputs
for i=1:nd
    fidi.rels(i).name  = dates(i).filename;   
end
for i=1:length(perms)
    fidi.perms(i).name = perms(i).filename;   
end
for i=1:length(c0s)
    fidi.c0(i).name    = c0s(i).filename;
end
%outputs
fido.resn0.name  = [rdir 'resn0'];
fido.shift.name  = [rdir 'shift'];
for i=1:nr
    fido.mag0(i).name  = [rdir rdates{i} '.mag0'];
    fido.time0(i).name = [rdir rdates{i} '.time0'];
    fido.magl(i).name  = [rdir rdates{i} '.maglow'];
    fido.magh(i).name  = [rdir rdates{i} '.maghigh'];
    fido.te(i).name    = [rdir rdates{i} '.timeerr'];
    fido.fwd5(i).name  = [rdir rdates{i} '.5day'];
    fido.status(i).name  = [rdir rdates{i} '.status'];
end
