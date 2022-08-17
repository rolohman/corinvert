function [dates,perms,c0,mag0,time0,output]=plot_all(x,y,pol,iscedir,fitdir,rlooks,alooks,dataflag,geoflag,plotflag,slcflag)
%iscedir - directory where isce was run (e.g., T130)
%fitDir - directory with the rain date fits. Empty to not plot. Should be
%on same grid as the first iscedir, geocoded.
%dataflag: 1: dir is isce, one track per directory. 2: iscedir and fitdir are resampled onto same grid (usually in overlap), multiple tracks.
%geoflag:  pixel given in: 1, geocoded pixels, 2, radar pixels, 3, latlon. for dataflag=2, geoflag = 3;
%plotflag: plot=1, noplot=0;
%slcflag:  save slc values = 1;

home=pwd;
relDirsuff=['/geo' pol];

if(or(rlooks>1,alooks>1))
    suff=['.' num2str(alooks) 'alks_' num2str(rlooks) 'rlks'];
else
    suff='';
end
if(~iscell(iscedir))
    iscedir={iscedir};
end

if(dataflag==2)
    geoflag=3;
    disp('x,y must be given in lon,lat for geoflag');
end

for i=1:length(iscedir)
    colf            = [iscedir{i} '/geo' pol '/cols' suff '.geo'];
    rowf            = [iscedir{i} '/geo' pol '/rows' suff '.geo'];
    
    switch geoflag
        case 1 %x, y row/col in first geocoded file, use same lat/lon for others, if present.
            
            if(i==1)
                xg(i)=x;
                yg(i)=y;
                latlonflag      = 2;
                [xr(i),yr(i),lon,lat] = LatLonRowCol(x,y,colf,rowf,latlonflag); %xr,yr in pixels, downlooked radar coords
                disp(['lon: ' num2str(lon) ' lat:' num2str(lat) ' col: ' num2str(round(xr(i))) ' row: ' num2str(round(yr(i)))])
                if(isnan(xr(i)))
                    disp('Must pick point in geocoded file that corresponds to data - you appear to have chosen a point outside extent of the image')
                    return
                end
            else
                latlonflag      = 1;
                [xr(i),yr(i),xg(i),yg(i)] = LatLonRowCol(lon,lat,colf,rowf,latlonflag); %xr,yr in pixels, downlooked radar coords
                disp(['lon: ' num2str(lon) ' lat:' num2str(lat) ' col: ' num2str(round(xr(i))) ' row: ' num2str(round(yr(i)))])
                if(isnan(xr(i)))
                    disp(['point chosen is outside of extent of track ' iscedir{i}])
                end
            end
            if(isfinite(xr(i)))
                chdir(iscedir{i})
                output  = plot_slc_rel(round(xr(i)),round(yr(i)),plotflag,slcflag);
                
                chdir(home);
            end
            
        case 2 %radar coords
            disp('havent checked this one yet.')
            
        case 3 %lat/lon
            if(dataflag==1)
                lon=x;
                lat=y;
                
                latlonflag      = 1;
                [xr(i),yr(i),xg(i),yg(i)] = LatLonRowCol(lon,lat,colf,rowf,latlonflag); %xr,yr in pixels, downlooked radar coords
                disp(['lon: ' num2str(lon) ' lat:' num2str(lat) ' col: ' num2str(round(xr(i))) ' row: ' num2str(round(yr(i)))])
                if(isnan(xr(i)))
                    disp(['point chosen is outside of extent of track ' iscedir{i}])
                end
                
                if(isfinite(xr(i)))
                    chdir(iscedir{i})
                    tmpout        = plot_slc_rel(round(xr(i)),round(yr(i)),plotflag,slcflag);
                    output(i,1:length(tmpout)) = tmpout;
                    chdir(home);
                end
            else
                files=dir([iscedir{i} '/rel*cor.geo']);
                file=[iscedir{i} '/' files(1).name];
                if(exist(file))
                    [~,b]=system(['grep rasterXSize ' file '.vrt']);
                    tmp=regexp(b,'rasterXSize="(\d+)" rasterYSize="(\d+)"','tokens');
                    if(length(tmp)==1)
                        nx=str2num(tmp{1}{1});
                        ny=str2num(tmp{1}{2});
                    else
                        disp('nx not found')
                    end
                    [~,b]=system(['grep GeoTransform ' file '.vrt']);
                    tmp=regexp(b,'[,<>]','split');
                    if(length(tmp)>1)
                        lon1=str2num(tmp{3});
                        dlon=str2num(tmp{4});
                        lat1=str2num(tmp{6});
                        dlat=str2num(tmp{8});
                        xg(i)  = floor((x-lon1)/dlon);
                        yg(i)  = floor((y-lat1)/dlat);
                    else
                        disp('GeoTransform not found')
                    end
                else
                    disp([file ' doesn''t exist']);
                end
                
            end
    end
end
xg
yg
if(fitdir)
    j=1; %assuming fitdir is on same grid as first iscedir for now.
    
    if(dataflag==1)
        relDir=[iscedir{j} '/geo' pol];
    elseif(dataflag==2)
        relDir = iscedir{j};
        rlooks=1;
        alooks=1;
    end
    %get data, rain dates, filehandles, sizes
    [dates,perms,c0s,dnr,fidi,fido,nxg,nyg]=pick_files_expfun(relDir,fitdir,rlooks,alooks,pol);
    dn    = [dates.dn];
    nd    = length(dates);
    nr    = length(dnr);
    nc    = length(c0s);
    
    if(xg(j)<1 || xg(j)>nxg || yg(j)<1 || yg(j)>nyg)
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
        fseek(fidi.c0(i).fid,(nxg*(yg(j)-1)+xg(j)-1)*4,-1);
    end
    if(fidi.shift.fid>0)
        fseek(fidi.shift.fid,(nxg*(yg(j)-1)+xg(j)-1)*4,-1);
    end
    for i=1:nr
        fseek(fidi.mag0(i).fid,(nxg*(yg(j)-1)+xg(j)-1)*4,-1);
        fseek(fidi.magl(i).fid,(nxg*(yg(j)-1)+xg(j)-1)*4,-1);
        fseek(fidi.magh(i).fid,(nxg*(yg(j)-1)+xg(j)-1)*4,-1);
        fseek(fidi.time0(i).fid,(nxg*(yg(j)-1)+xg(j)-1)*4,-1);
        fseek(fidi.te(i).fid,(nxg*(yg(j)-1)+xg(j)-1)*4,-1);
    end
    for i=1:nd
        fseek(fidi.rels(i).fid,(nxg*(yg(j)-1)+xg(j)-1)*4,-1);
    end
    for i=1:length(perms)
        fseek(fidi.perms(i).fid,(nxg*(yg(j)-1)+xg(j)-1)*4,-1);
    end
    
    for i=1:nc
        c0(i)    = fread(fidi.c0(i).fid,1,'real*4');
    end
    c0=median(c0,'omitnan');
    if(fidi.shift.fid>0)
        shift = fread(fidi.shift.fid,1,'real*4');
    else
        shift=0;
    end
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
end
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

