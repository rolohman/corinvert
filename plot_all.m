function [dates,perms,c0,mag0,time0,output]=plot_all(x,y,pol,relDir,rdir,iscedir,rlooks,alooks,geoflag,plotflag,slcflag)
%geoflag:  pixel given in: 1, geocoded pixels, 2, radar pixels, 3, latlon
%plotflag: plot=1, noplot=0;
%slcflag:  save slc values = 1;
%for now gflag=1;
home=pwd;

if(geoflag==1)
    %input coords are from geocoded file.
    xg              = x;
    yg              = y;
    if(or(rlooks>1,alooks>1))
        suff=['.' num2str(alooks) 'alks_' num2str(rlooks) 'rlks'];
    else
        suff='';
    end
    
    if(length(iscedir)==1) %one track only
        latlonflag      = 2;
        colf            = [relDir '/cols' suff '.geo'];
        rowf            = [relDir '/rows' suff '.geo'];
        [xr,yr,lon,lat] = LatLonRowCol(x,y,colf,rowf,latlonflag); %xr,yr in pixels, downlooked radar coords
        disp(['lon: ' num2str(lon) ' lat:' num2str(lat)])
        chdir(iscedir)
        [output]        = plot_slc_rel(round(xr),round(yr),plotflag,slcflag);
        chdir(home);
    elseif(length(iscedir)>1) %multiple tracks, resampled to same grid in dir, iscedir-> orig geo dirs
        latlonflag      = 2;
        tmp             = dir([relDir '/*.geo']);
        tmpf            = [tmp(1).folder '/' tmp(1).name];
        [~,~,lon,lat]   = LatLonRowCol(x,y,tmpf,tmpf,latlonflag); %xr,yr in pixels, downlooked radar coords
      
        
        for i=1:length(iscedir)
         
            [xr,yr,lon,lat] = LatLonRowCol(x,y,colf,rowf,latlonflag); %xr,yr in pixels, downlooked radar coords
         
                colf            = [iscedir{i} '/geo' pol '/cols' suff '.geo'];
            rowf            = [iscedir{i} '/geo' pol '/rows' suff '.geo'];
        
            disp(['lon: ' num2str(lon) ' lat:' num2str(lat)])
            chdir(iscedir)
            [output]        = plot_slc_rel(round(xr),round(yr),plotflag,slcflag);
            chdir(home);
        end
        
    end
else
    disp('gflag option not done yet for gflag ne 1')
    return
end

%get data, rain dates, filehandles, sizes
[dates,perms,c0s,dnr,fidi,fido,nxg,nyg]=pick_files_expfun(relDir,rdir,rlooks,alooks,pol);
dn    = [dates.dn];
nd    = length(dates);
nr    = length(dnr);
nc    = length(c0s);

if(xg<1 || xg>nxg || yg<1 || yg>nyg)
    disp('point out of range')
    return;
end

%open all geocoded files
names=fieldnames(fido); %we are opening these all as read/input
for i=1:length(names)
    fidi.(names{i})=fido.(names{i});
end
[fidi,~,~]=open_files(fidi,[],nxg,nyg); %don't want to overwrite outfiles!

%get data at points in geocoded files
for i=1:nc
    fseek(fidi.c0(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
end
fseek(fidi.shift.fid,(nxg*(yg-1)+xg-1)*4,-1);
for i=1:nr
    fseek(fidi.mag0(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
    fseek(fidi.magl(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
    fseek(fidi.magh(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
    fseek(fidi.time0(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
    fseek(fidi.te(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
end
for i=1:nd
    fseek(fidi.rels(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
end
for i=1:length(perms)
    fseek(fidi.perms(i).fid,(nxg*(yg-1)+xg-1)*4,-1);
end

for i=1:nc
    c0(i)    = fread(fidi.c0(i).fid,1,'real*4');
end
c0=median(c0,'omitnan');
shift = fread(fidi.shift.fid,1,'real*4');
for i=1:length(dnr)
    mag0(i)    = -log(fread(fidi.mag0(i).fid,1,'real*4'));
    magl(i)    = -log(fread(fidi.magl(i).fid,1,'real*4'));
    magh(i)    = -log(fread(fidi.magh(i).fid,1,'real*4'));
    time0(i)   = fread(fidi.time0(i).fid,1,'real*4');
    timeerr(i) = fread(fidi.te(i).fid,1,'real*4');
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
goodmod = find(isfinite(mag0+time0));
timemat=zeros(length(dn2),length(goodmod)); 
for i=1:length(goodmod)
    timemat(:,i) = dn2-dnr(goodmod(i));
end
timemat(timemat<0)=0;

mod0      = [time0(goodmod) mag0(goodmod) shift];
synth     = exp(expfun_all(mod0,timemat,zeros(size(dn2')))); %negative since fun outputs residual, y-synth.
for i=1:nd
    dates(i).synth=synth(i);
end
synthnan    = [synth;nan*dnr];
[dnnan,sid] = sort([dn2 dnr']); %add points to break for nans
synthnan    = synthnan(sid);

%plotting

if(plotflag)
    figure('Name',[num2str(xg) ' ' num2str(yg)])
    
    plot([dnr dnr]',exp(-[magl;magh]),'-','color',[.4 .4 .4],'linewidth',3)
    hold on
    plot([dates.dn],[dates.rel],'.','color','b','markersize',12)
    
    for j=1:length(perms)
        plot([perms(j).d1 perms(j).d2],[perms(j).perm perms(j).perm],'-','color','r','linewidth',2)
    end
    
    plot(dnnan,synthnan,'k.-') %nans and sorting make a discontinuous line, broken at dates.
    
    axis tight,ax=axis;axis([ax(1:2) 0 ax(4)]);ax=axis;
    datetick('x','yymmm')
    for i=1:length(dnr)
        plot([dnr(i) dnr(i)],ax(3:4),'m--');
    end
    %plot(ax(1:2),[shift shift],'g:')
    grid on
    axis(ax)

end

