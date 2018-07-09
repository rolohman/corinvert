function [dates,perms]=plot_fit_latlon(xpt,ypt)
dirs={'156_674','47_1096','54_673'};
dirs={'156_674','54_673'};
nd=0;
dates=[];
perms=[];
for i=1:length(dirs)
    tmp=dir([dirs{i} '/geotiffs/geofiles/rel*.cor']);
    
    for j=1:length(tmp)
        t=tmp(j).name;
        dates(end+1).name=[dirs{i} '/geotiffs/geofiles/' t];
        dates(end).dn=datenum(t(5:12),'yyyymmdd');
        dates(end).t=i;
    end
    tmp=dir([dirs{i} '/geotiffs/geofiles/perm*cor']);
    
    for j=1:length(tmp)
        t=tmp(j).name;
        perms(end+1).name=[dirs{i} '/geotiffs/geofiles/' t];
        perms(end).d1=datenum(t(6:13),'yyyymmdd');
        perms(end).d2=datenum(t(15:22),'yyyymmdd');
        perms(end).t=i;
        perms(end).d=mean([perms(end).d1 perms(end).d2]);
    end
end

d=[dates.t];
nd=length(dates);
dn    = [dates.dn];

nx_geo    = 7560;
ny_geo    = 6228;
x1=-70.6;
dx=0.000277777777777777;
y1=-24.0;
dy=-0.000277777777777778;

%
disp(x1+xpt*dx)
disp(y1+ypt*dy)

if(xpt<1 || xpt>nx_geo || ypt<1 || ypt>ny_geo)
    disp('point out of range')
    return;
end


rdir = 'results_TS/';
rdate  = {'20150325','20150807','20160625','20170310','20170513','20170606'};

dnr    = datenum(rdate,'yyyymmdd');
rdate  = rdate(dnr<max(dn));
dnr    = dnr(dnr<max(dn));


for i=1:nd
    fid=fopen(dates(i).name,'r');
    fseek(fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
    dates(i).rel=fread(fid,1,'real*4');
    fclose(fid);
    if(dates(i).rel==-9999)
        dates(i).rel=NaN;
    end
end

for i=1:length(perms)
    fid=fopen(perms(i).name,'r');
    fseek(fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
    perms(i).perm=fread(fid,1,'real*4');
    fclose(fid);
    if(perms(i).perm==-9999)
        perms(i).perm=NaN;
    end
end

[jnk,sid]=sort([dates.dn]);


tc={'r','b','g'};
figure
hold on
for i=1:3
    id1=find([dates.t]==i);
    plot([dates(id1).dn],[dates(id1).rel],'.','color',tc{i},'markersize',12)
    id1=find([perms.t]==i);

     for j=1:length(id1)
         plot([perms(id1(j)).d1 perms(id1(j)).d2],[perms(id1(j)).perm perms(id1(j)).perm],'-','color',tc{i})
     end
end



axis tight
ax=axis;
datetick('x','yymmm')
for i=1:length(dnr)
    plot([dnr(i) dnr(i)],ax(3:4),'m--');
end
grid on
axis(ax)


%

