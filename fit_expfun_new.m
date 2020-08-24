addpath('~/matlab/DERIVESTsuite');
parpool(10)
home=pwd;
dirs={'T130'};
%dirs={'T28'};
%dirs={'T101'};
pol='_VV';
corcutoff=0.01;
options = optimset('Display','off','TolFun',1e-4);

nt=length(dirs);
nd=0;
dates=[];
for i=1:length(dirs)
    tmp=dir([dirs{i} '/geo' pol '/rel*_4r_4a.cor.geo']);
    for j=1:length(tmp)
        t=tmp(j).name;
        vrt=[tmp(j).folder '/' t '.vrt'];
        dates(end+1).name=[tmp(j).folder '/' t];
        dates(end).dn=datenum(t(5:12),'yyyymmdd');
    end
end

%get nx/ny

if(exist(vrt))
    [a,b]=system(['grep rasterXSize ' vrt]);
    tmp=regexp(b,'rasterXSize="(\d+)" rasterYSize="(\d+)">','tokens');
    if(length(tmp)==1)
        nx=str2num(tmp{1}{1});
        ny=str2num(tmp{1}{2});
        [nx ny]
    end
else
    disp('no nx/ny info, need vrt file')
    return
end

[jnk,sortid]=sort([dates.dn]);
dates = dates(sortid);
dn    = [dates.dn];
nd    = length(dates);


rdir = ['results_TS' pol '/'];
rdate  = {'20150325','20150807','20160625','20170310','20170513','20170606'};
rdate = {'20170325','20170506','20170717','20171009','20171220','20180302','20180524','20181011'};
if(~exist(rdir,'dir'))
    mkdir(rdir)
end
dnr    = datenum(rdate,'yyyymmdd');
goodr   = and(dnr<max(dn),dnr>min(dn)); %only use rain during time series
rdate  = rdate(goodr);
dnr    = dnr(goodr);
nr     = length(dnr);
  
timemat=zeros(nd,nr); %used later in inversion.
for i=1:nr
    timemat(:,i)=dn-dnr(i);
    befid(i)=find(dn<dnr(i),1,'last');
    afid(i)=find(dn>dnr(i),1,'first');
end
dt = dn(afid)-dnr'; %time since event, used in res.
timemat(timemat<0)=0;
LB     = dt+3; %lower bound on time - at least 1/2 time after first date.
UB     = 100*ones(1,nr); %upper bound on time
tmod0  = 12*ones(1,nr); %starting time model

%%  Open files
for i=1:nd
    fidr(i)=fopen(dates(i).name,'r');
end
for i=1:length(dirs)
    tmp=[dirs{i} '/geo' pol '/c0_4r_4a.cor.geo'];
    fid0(i)=fopen(tmp,'r');
end
%find last line written
testfile=[rdir rdate{1} '.mag0'];
if(exist(testfile,'file'))
    fid=fopen(testfile,'r');
    for j=1:ny
        [tmp,count]=fread(fid,nx,'real*4');
        if(count<nx)
            break
        end
    end
    fclose(fid);
    
    online=j-1
    if(online<1)
        disp([ testfile  ' empty?'])
        return
    end
    for i=1:nr
        fidmag(i)  = fopen([rdir rdate{i} '.mag0'],'a');
        fidt(i)    = fopen([rdir rdate{i} '.time0'],'a');
        fidmagl(i) = fopen([rdir rdate{i} '.maglow'],'a');
        fidmagh(i) = fopen([rdir rdate{i} '.maghigh'],'a');
        fidte(i)   = fopen([rdir rdate{i} '.timeerr'],'a');
        fidmod5(i) = fopen([rdir rdate{i} '.5day'],'a');
    end
    
    fidres  = fopen([rdir 'resn0'],'a');
    fids    = fopen([rdir 'shift'],'a');
    
    for i=1:nd
        fseek(fidr(i),online*nx*4,-1);
    end
    for i=1:length(dirs)
        fseek(fid0(i),online*nx*4,-1);
    end
else
    online=0;
    for i=1:nr
        fidmag(i)  = fopen([rdir rdate{i} '.mag0'],'w');
        fidt(i)    = fopen([rdir rdate{i} '.time0'],'w');
        fidmagl(i) = fopen([rdir rdate{i} '.maglow'],'a');
        fidmagh(i) = fopen([rdir rdate{i} '.maghigh'],'a');
        fidte(i)   = fopen([rdir rdate{i} '.timeerr'],'a');
        fidmod5(i) = fopen([rdir rdate{i} '.5day'],'a');
    end
    fidres  = fopen([rdir 'resn0'],'w');
    fids    = fopen([rdir 'shift'],'w');
end


%% start
for j=online+1:ny
    j
    dat=zeros(nd,nx);
    for i=1:nd
        [tmp,count]=fread(fidr(i),nx,'real*4');
        if(count>0)
            dat(i,1:count)=tmp;
        end
    end
    bad      = or(dat==-9999,or(isinf(dat),dat==0));
    dat(bad) = NaN;
    count    = sum(isfinite(dat));
    sig      = -0.1*dat.^2+0.01*dat+.09; %based on fit to numerical 60-look cor estimate, cor>0.3.
    sig      = max(sig,corcutoff); %can't really determine cor to within cutoff
    dhlog    = -log(dat-sig);
    dllog    = -log(dat+sig);
    dsig     = (dhlog-dllog)/2; % approx std. dev of log values.
    dsig(dat<sig)=dllog(dat<sig); %removes imaginary values
    
    c0s=nan(length(dirs),nx);
    for i=1:length(dirs)
        [tmp,count2]=fread(fid0(i),nx,'real*4');
        if(count2>0)
            c0s(i,1:count2)=tmp;
        end
    end
    bad      = or(c0s==-9999,or(isinf(c0s),c0s==0));
    c0s(bad) = NaN;
    c0s      = median(c0s,1,'omitnan');
    
    goodid   = find(and(count>20,c0s>0.3));
    
    mags    = zeros(nr,length(goodid));
    times   = nan*mags;
    maglow  = mags;
    maghig  = mags;
    timerr  = mags;
    mag5    = mags;
    shifts  = nan(1,length(goodid));
    allres  = shifts;
    
    parfor i=1:length(goodid)
        
        logd      = -log(dat(:,goodid(i)));
        deld      = logd(afid)-logd(befid);
        
        goodd     = isfinite(logd);
        weights   = diag(dsig(goodd,goodid(i)).^2); %data covariance matrix
        y         = logd(goodd);
        magi      = zeros(nr,1);
        goodr     = deld>corcutoff;
        ng        = sum(goodr);
        notdone   = ng>0;
        tmptime   = NaN(nr,1);
        tmpmag    = NaN(nr,1);
        tstd      = [];
        mstd      = [];
        res       = [];
         shift     = 0;
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
                   
                end
            end
            if(or(ng==0,toobad==nr))
                tmpmag(:)           = NaN;
                tmptime(:)          = NaN;
                tstd                = [];
                mstd                = [];
                notdone             = 0;
                shift               = 0;    
            end
        end
        
        tmpmag5=tmpmag.*exp(-5./tmptime);
   
        times(:,i) = tmptime;
        mags(:,i)  = tmpmag;
        mag5(:,i)  = tmpmag5;
        
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

    end
    
    tmp=nan(1,nx);
    for i=1:nr
        tmp(goodid) = exp(-mags(i,:));
        fwrite(fidmag(i),tmp,'real*4');
        tmp(goodid) = times(i,:);
        fwrite(fidt(i),tmp,'real*4');
        tmp(goodid) = exp(-maglow(i,:));
        fwrite(fidmagl(i),tmp,'real*4');
        tmp(goodid) = exp(-maghig(i,:));
        fwrite(fidmagh(i),tmp,'real*4');
        tmp(goodid) = timerr(i,:);
        fwrite(fidte(i),tmp,'real*4');
        tmp(goodid) = exp(-mag5(i,:));
        fwrite(fidmod5(i),tmp,'real*4');
    end
    tmp(goodid) = allres;
    fwrite(fidres,tmp,'real*4');
    tmp(goodid) = shifts;
    fwrite(fids,tmp,'real*4');

end
fclose('all');
