addpath('~/matlab/DERIVESTsuite');
parpool(10)
home=pwd;
relDir='T130';
%relDir='T28';
%relDir='T101';
pol='_VV';
corcutoff      = 0.01; %diff to distinguish from 1
lowcorcutoff   = 0.3; %lowest background cor (c0)
lowcountcutoff = 20;  %need at least n dates

options = optimset('Display','off','TolFun',1e-4);

dates = [];
tmp=dir([relDir '/geo' pol '/rel*_4r_4a.cor.geo']);
for j=1:length(tmp)
    t=tmp(j).name;
    dates(end+1).name   = t(5:12);
    dates(end).filename = [tmp(j).folder '/' t];
    dates(end).dn       = datenum(t(5:12),'yyyymmdd');
    dates(end).origdir  = [tmp(j).folder];
end

%get nx/ny from last file opened
vrt=[dates(end).filename '.vrt'];
if(exist(vrt,'file'))
    [a,b] = system(['grep rasterXSize ' vrt]);
    tmp   = regexp(b,'rasterXSize="(\d+)" rasterYSize="(\d+)">','tokens');
    if(length(tmp)==1)
        nx=str2num(tmp{1}{1}); ny=str2num(tmp{1}{2});
        [nx ny]
    end
else
    disp('no nx/ny info, need vrt file'),return
end

[jnk,sortid]=sort([dates.dn]);
dates = dates(sortid);
dn    = [dates.dn];
nd    = length(dates);

rdir = ['results_TS' pol '/'];
rdate=cellstr(num2str(load('dates.list')));
if(~exist(rdir,'dir'))
    mkdir(rdir)
end

dnr    = datenum(rdate,'yyyymmdd');
goodr  = and(dnr<max(dn),dnr>min(dn)); %only use rain during time series
rdate  = rdate(goodr);
dnr    = dnr(goodr);
nr     = length(dnr);
  
timemat=zeros(nd,nr); %used later in inversion.
for i=1:nr
    timemat(:,i) = dn-dnr(i);
    befid(i)     = find(dn<dnr(i),1,'last');
    afid(i)      = find(dn>dnr(i),1,'first');
end
timemat(timemat<0)=0;
dt     = dn(afid)-dnr'; %time since event, used in res.
LB     = dt+3; %lower bound on time - at least 1/2 time after first date.
UB     = 100*ones(1,nr); %upper bound on time
tmod0  = 12*ones(1,nr); %starting time model

%open all filesif(type==1)
clear fidi fido
%inputs
for i=1:nd
    fidi.rels(i).name=dates(i).filename;   
end
fidi.c0.name     = [relDir '/geo' pol '/c0_4r_4a.cor.geo'];
%outputs
fido.resn0.name  = [rdir 'resn0'];
fido.shift.name  = [rdir 'shift'];
for i=1:nr
    fido.mag0(i).name  = [rdir rdate{i} '.mag0'];
    fido.time0(i).name = [rdir rdate{i} '.time0'];
    fido.magl(i).name  = [rdir rdate{i} '.maglow'];
    fido.magh(i).name  = [rdir rdate{i} '.maghigh'];
    fido.te(i).name    = [rdir rdate{i} '.timeerr'];
    fido.fwd5(i).name  = [rdir rdate{i} '.5day'];
    fido.status(i).name  = [rdir rdate{i} '.status'];
end

[fidi,fido,online]=open_files(fidi,fido,nx,ny);


%% start
for j=online+1:ny
    j
    dat=zeros(nd,nx);
    for i=1:nd
        [tmp,count1]=fread(fidi.rels(i).fid,nx,'real*4');
        if(count1>0)
            dat(i,1:count1)=tmp;
        end
    end
    bad      = or(dat==-9999,or(isinf(dat),dat==0));
    dat(bad) = NaN;
    count    = sum(isfinite(dat)); %how many dates have observations
    
    %calculate expected error on corr measurements (larger for low cor)
    sig      = -0.1*dat.^2+0.01*dat+.09; %based on fit to numerical 60-look cor estimate, cor>0.3.
    sig      = max(sig,corcutoff); %can't really determine cor to within cutoff
    
    %translate errors to log space, high and low, diff.
    dhlog    = -log(dat-sig); dllog    = -log(dat+sig);
    dsig     = (dhlog-dllog)/2; % approx std. dev of log values.
    dsig(dat<sig)=dllog(dat<sig); %removes imaginary values
    
    %load "background" c0 values - won't use pixels with really low values.
    c0s=nan(1,nx);
    [tmp,count2]=fread(fidi.c0.fid,nx,'real*4');
    if(count2>0)
        c0s(1:count2)=tmp;
    end 
    c0s(or(c0s==-9999,or(isinf(c0s),c0s==0))) = NaN;   
    
    goodid   = find(and(count>lowcountcutoff,c0s>lowcorcutoff));
    
    %initialize arrays (necessary for parfor)
    mags    = nan(nr,length(goodid));times = mags; maglow = mags; maghig = mags; timerr = mags; fwd5 = mags; status=mags;
    shifts  = nan(1,length(goodid)); allres = shifts;
    
    parfor i=1:length(goodid)
        
        logd      = -log(dat(:,goodid(i)));
        deld      = logd(afid)-logd(befid);
        
        goodd     = isfinite(logd);
        weights   = diag(dsig(goodd,goodid(i)).^2); %data covariance matrix
        y         = logd(goodd);
        magi      = zeros(nr,1);
        goodr     = deld>corcutoff;
        ng        = sum(goodr);
        notdone   = ng>0; %start loop, as long as there are some "good" values
       
        
        %%%initialize values
        tmptime = NaN(nr,1);tmpstat=ones(nr,1);tmpmag = tmptime;tstd = [];mstd = [];res = [];shift = 0; %mark tstd with 20 to see where ng cutoff occurrs
        if(ng==0)
            tmpstat(~goodr)=2; %no rain dates with pos cutoff
        end
        while(notdone)
            x                   = timemat(goodd,goodr);
            [mod1,~,~,~,~,~,J0] = lsqnonlin('expfun',tmod0(goodr),LB(goodr),UB(goodr),options,x,y);

            jt                  = sum(abs(full(J0)));           
            [res,magst,shift,G] = expfun(mod1,x,y);
            
            magi(goodr)         = magst;
            magi(~goodr)        = 0;
            
            syntmp              = magst'.*exp(-dt(goodr)./mod1); %pred cor at dates after event
            
            toobad              = or(jt==0,syntmp<=corcutoff);
            goodr(goodr)        = ~toobad;
            toobad              = sum(toobad);
            ng                  = sum(goodr);
            if(toobad==0)
        
                mod0                = [mod1 magst' shift];
                J                   = jacobianest('expfun_all',mod0,x,y);
                J2                  = J'*J;
                
                allJ(i)             = rcond(J2);
                J3                  = inv(J2);
                modcov              = J3*J'*weights*J*J3;
                modstd              = sqrt(diag(modcov));
                
                tstd                = modstd(1:ng);
                mstd                = modstd(ng+[1:ng]);
                
                badt                = tstd>100;
                if(sum(badt))
                    notdone         = 1;
                    gid=find(goodr);
                    goodr(gid(badt))=false;
                    ng=sum(goodr);
                     
                else
                    notdone             = 0;
                    tmptime(goodr)      = mod1;
                    tmpmag(goodr)       = magst';
                    tmpstat(goodr)      = 3; %worked!
                end
            end
            if(or(ng==0,toobad==nr))
                tmpmag(:)           = NaN;
                tmptime(:)          = NaN;
                tstd                = [];
                mstd                = [];
                notdone             = 0;
                shift               = 0;    
                tmpstat(:)          = 4; %all bad?
            end
        end
        
        tmpmag5=tmpmag.*exp(-5./tmptime);
   
        %below is bookkeeping related to the use of parfor loops.
        times(:,i) = tmptime;
        mags(:,i)  = tmpmag;
        fwd5(:,i)  = tmpmag5;
        
        tmp=nan(nr,1);
        tmp(goodr)      = tstd;
        timerr(:,i)     = tmp;
        magstd          = mstd;
        tmp             = nan(nr,1);
        tmp(goodr)      = magi(goodr)+magstd;
        maglow(:,i)     = tmp;
        tmp             = nan(nr,1);
        tmp(goodr)      = magi(goodr)-magstd;
        maghig(:,i)     = tmp;
        allres(i)       = std(res);
        shifts(i)       = shift;
        status(:,i)     = tmpstat;
    end
    
    tmp=nan(1,nx);
    for i=1:nr
        tmp(goodid) = exp(-mags(i,:));
        fwrite(fido.mag0(i).fid,tmp,'real*4');
        tmp(goodid) = times(i,:);
        fwrite(fido.time0(i).fid,tmp,'real*4');
        tmp(goodid) = exp(-maglow(i,:));
        fwrite(fido.magl(i).fid,tmp,'real*4');
        tmp(goodid) = exp(-maghig(i,:));
        fwrite(fido.magh(i).fid,tmp,'real*4');
        tmp(goodid) = timerr(i,:);
        fwrite(fido.te(i).fid,tmp,'real*4');
        tmp(goodid) = exp(-fwd5(i,:));
        fwrite(fido.fwd5(i).fid,tmp,'real*4');
        tmp(goodid) = status(i,:);
        fwrite(fido.status(i).fid,tmp,'integer*1');
    end
    tmp(goodid) = allres;
    fwrite(fido.resn0.fid,tmp,'real*4');
    tmp(goodid) = shifts;
    fwrite(fido.shift.fid,tmp,'real*4');

end
fclose('all');
