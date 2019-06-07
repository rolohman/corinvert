function [dates,perms,c0]=plot_fit_latlon(xpt,ypt,pol,pflag)
dirs={'T54','T156','T47'};
%dirs=dirs(1:2);
nd=0;
dates=[];
perms=[];
c0=[];
nx_geo    = 2340;
ny_geo    = 2574;
for i=1:length(dirs)
    tmp=dir([dirs{i} '/geo' pol '/rel*_4r_4a.cor.geo']);
    
    for j=1:length(tmp)
        t=tmp(j).name;
        dates(end+1).name=[tmp(j).folder '/' t];
        dates(end).dn=datenum(t(5:12),'yyyymmdd');
        dates(end).t=i;
    end
    tmp=dir([dirs{i}  '/geo' pol '/perm*_4r_4a.cor.geo']);
    
    for j=1:length(tmp)
        t=tmp(j).name;
        perms(end+1).name=[tmp(j).folder '/' t];
        perms(end).d1=datenum(t(6:13),'yyyymmdd');
        perms(end).d2=datenum(t(15:22),'yyyymmdd');
        perms(end).t=i;
        perms(end).d=mean([perms(end).d1 perms(end).d2]);
    end
    file=[dirs{i}  '/geo' pol '/c0_4r_4a.cor.geo'];
    fid=fopen(file,'r');
    fseek(fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
    c0(end+1)=  fread(fid,1,'real*4');
    fclose(fid);
end

d=[dates.t];
nd=length(dates);
dn    = [dates.dn];


% x1=-70.6;
% dx=0.000277777777777777;
% y1=-24.0;
% dy=-0.000277777777777778;
% 
% %
% disp(x1+xpt*dx)
% disp(y1+ypt*dy)

if(xpt<1 || xpt>nx_geo || ypt<1 || ypt>ny_geo)
    disp('point out of range')
    return;
end


rdir = ['results_TS' pol '/'];
rdate  = {'20150325','20150807','20160625','20170310','20170513','20170606'};

dnr    = datenum(rdate,'yyyymmdd');
rdate  = rdate(dnr<max(dn));
dnr    = dnr(dnr<max(dn));

dn2    = [dn dnr']; %add points at times of rain for more complete plotting
synth  = 0*dn2;
for i=1:length(rdate)
    file=[rdir rdate{i} '.mag0'];
    fid=fopen(file,'r');
    fseek(fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
    mag=-log(fread(fid,1,'real*4'));
    file=[rdir rdate{i} '.maglow'];
    fid=fopen(file,'r');
    fseek(fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
    magl(i)=fread(fid,1,'real*4');
    file=[rdir rdate{i} '.maghigh'];
    fid=fopen(file,'r');
    fseek(fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
    magh(i)=fread(fid,1,'real*4');
    file=[rdir rdate{i} '.time0'];
    fid=fopen(file,'r');
    fseek(fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
    time=fread(fid,1,'real*4');
    
    id=find(dn2>=dnr(i));
    if(and(isfinite(mag),isfinite(time)))
        synth(id)=synth(id)+mag*exp(-(dn2(id)-dnr(i))/time);
    end
end
synth=exp(-synth);

for i=1:nd
    dates(i).synth=synth(i);
    fid=fopen(dates(i).name,'r');
    fseek(fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
    dates(i).rel=fread(fid,1,'real*4');
    fclose(fid);
    if(or(dates(i).rel==-9999,dates(i).rel==0))
        dates(i).rel=NaN;
    end
end

for i=1:length(perms)
    fid=fopen(perms(i).name,'r');
    fseek(fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
    perms(i).perm=fread(fid,1,'real*4');
    fclose(fid);
    if(or(perms(i).perm==-9999,perms(i).perm==0))
        perms(i).perm=NaN;
    end
end

dn2=[dn2 dnr'-0.001]; %add points to break for nans
synth2=[synth nan*dnr'];
[jnk,sid]=sort(dn2);

if(plag)
tc={'r','[0.1 0 0.9]','[0 0.6 0]'};
figure('Name',[num2str(xpt) ' ' num2str(ypt)])

plot([dnr dnr]',[magl;magh],'-','color',[.4 .4 .4],'linewidth',3)
hold on
for i=1:length(dirs)
    id1=find([dates.t]==i);
    plot([dates(id1).dn],[dates(id1).rel],'.','color',tc{i},'markersize',12)
    id1=find([perms.t]==i);

     for j=1:length(id1)
         plot([perms(id1(j)).d1 perms(id1(j)).d2],[perms(id1(j)).perm perms(id1(j)).perm],'-','color',tc{i},'linewidth',2)
     end
end

plot(dn2(sid),synth2(sid),'k.-') %nans and sorting make a discontinuous line, broken at dates.


axis tight
ax=axis;
axis([ax(1:2) 0 ax(4)]);
ax=axis;
datetick('x','yymmm')
for i=1:length(dnr)
    plot([dnr(i) dnr(i)],ax(3:4),'m--');
end
grid on
axis(ax)

end
%

