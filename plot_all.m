function [dates,perms,c0,mag0,time0]=plot_all(x,y,pol,relDir,gflag,pflag)
%gflag: pixel given in: 1, geocoded pixels, 2, radar pixels, 3, latlon
%for now gflag=1;
if(gflag==1)
    %input coords are from geocoded file.
    xg         = x;
    yg         = y;
    latlonflag = 2; 
    colf       = [relDir '/geo' pol '/cols.4alks_4rlks.cor.geo'];
    rowf       = [relDir '/geo' pol '/rows.ralks_4rlks.cor.geo'];
    [xr,yr,lon,lat] = LatLonRowCol(x,y,colf,rowf,latlonflag); %xr,yr in pixels, downlooked radar coords
else
    disp('gflag option not done yet for gflag ne 1')
    return
end

%get data, rain dates, filehandles, sizes
[dates,perms,rdates,dnr,fidi,fido,nxg,nyg]=pick_files_expfun(relDir,pol);
dn    = [dates.dn];
nd    = length(dates);
nr    = length(dnr);

if(xg<1 || xg>nxg || yg<1 || yg>nyg)
    disp('point out of range')
    return;
end

%open all geocoded files
names=fieldnames(fido); %we are opening these all as read/input
for i=1:length(names)
    fidi.(names{i})=fido.(names{i});
end
[fidi,~,online]=open_files(fidi,[],nxg,nyg); %don't want to overwrite outfiles!
disp(['processed expfit to ' num2str(online)])


%get data at points in geocoded files
fseek(fidi.c0.fid,(nxg*(yg-1)+xg-1)*4,-1);
fseek(fidi.shift.fid,(nxg*(yg-1)+xg-1)*4,-1);
for i=1:length(dnr)
    fseek(fidi.mag0(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
    fseek(fidi.magl(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
    fseek(fidi.magh(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
    fseek(fidi.time0(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
    fseek(fidi.timeerr(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
end
for i=1:nd
    fseek(fidi.rels(i).fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
end
for i=1:length(perms)
    fseek(fidi.perms(i).fid,(nx_geo*(ypt-1)+xpt-1)*4,-1);
end

c0    = fread(fidi.c0.fid,1,'real*4');
shift = fread(fidi.shift.fid,1,'real*4');
for i=1:length(dnr)
    mag0(i)    = -log(fread(fidi.mag0(i).fid,1,'real*4'));
    magl(i)    = -log(fread(fidi.magl(i).fid,1,'real*4'));
    magh(i)    = -log(fread(fidi.magh(i).fid,1,'real*4'));
    time0(i)   = fread(fidi.time0(i).fid,1,'real*4');
    timeerr(i) = fread(fidi.timeerr(i).fid,1,'real*4');
end
for i=1:nd
    dates(i).rel=fread(fidi.rels(i).fid,1,'real*4');
    if(or(dates(i).rel==-9999,dates(i).rel==0))
        dates(i).rel=NaN;
    end
end
for i=1:length(perms)
    perms(i).perm=fread(fidi.perms(i).fid,1,'real*4');
    if(or(perms(i).perm==-9999,perms(i).perm==0))
        perms(i).perm=NaN;
    end
end
fclose('all');

%%%synthetic
dn2    = [dn dnr'+0.01]; %add points at times of rain for more complete plotting
%Matrix of times since rain, used in forward model
timemat=zeros(length(dnr2),nr); 
for i=1:nr
    timemat(:,i) = dn2-dnr(i);
end
timemat(timemat<0)=0;

mod0      = [time0 mag0 shift];
synth     = exp(expfun_all(mod0,timemat,zeros(size(dn2)))); %negative since fun outputs residual, y-synth.
for i=1:nd
    dates(i).synth=synth(i);
end
dn2       = [dn2 dnr']; %add points to break for nans
synth2    = [synth nan*dnr'];
[dn3,sid] = sort(dn2);
synth2    = synth2(sid);

%plotting

if(pflag)
    figure('Name',[num2str(xg) ' ' num2str(yg)])
    
    plot([dnr dnr]',[magl;magh],'-','color',[.4 .4 .4],'linewidth',3)
    hold on
    plot([dates.dn],[dates.rel],'.','color','b','markersize',12)
    
    for j=1:length(perms)
        plot([perms(j).d1 perms(j).d2],[perms(j).perm perms(j).perm],'-','color','r','linewidth',2)
    end
    
    plot(dn3,synth3,'k.-') %nans and sorting make a discontinuous line, broken at dates.
    
    axis tight,ax=axis;axis([ax(1:2) 0 ax(4)]);ax=axis;
    datetick('x','yymmm')
    for i=1:length(dnr)
        plot([dnr(i) dnr(i)],ax(3:4),'m--');
    end
    grid on
    axis(ax)

end

