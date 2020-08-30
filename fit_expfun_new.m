addpath('~/matlab/DERIVESTsuite');
home=pwd;
%% check directories
if(~exist('pol','var'))
    disp('should define pol, using _VV');
pol='_VV';
else
    disp(['using polarity ' pol])
end
if(~exist('relDir','var'))
    disp('must define relDir, full path to geocoded relative coherence values')
    return
else
    disp(['using coherence in ' relDir])
end
if(~exist('rdir','var'))
    disp('must define output directory')
    return
else
    disp(['output in ' rdir])
    if(~exist(rdir,'dir'))
        mkdir(rdir)
    end
end
if(~exist('rlooks','var'))
    rlooks=4;
end
if(~exist('alooks','var'))
    alooks=4;
end

%% define cutoffs
corcutoff      = 0.01; %diff to distinguish from 1
lowcorcutoff   = 0.35; %lowest background cor (c0)
lowcountcutoff = 20;  %need at least n dates
maxT           = 100; %maxmum time for exponential fit, days
startT         = 12;  %time for each event, initialize, days
options        = optimset('Display','off','TolFun',1e-4);

%% get data, rain dates, filehandles, sizes
[dates,perms,c0s,rdates,dnr,fidi,fido,nx,ny]=pick_files_expfun(relDir,rdir,rlooks,alooks,pol);
dn    = [dates.dn];
nd    = length(dates);
nr    = length(dnr);
nc    = length(c0s);

%% Matrix of times since rain, used later in inversion.
timemat=zeros(nd,nr); 
for i=1:nr
    timemat(:,i) = dn-dnr(i);
    befid(i)     = find(dn<dnr(i),1,'last'); %first date before rain
    afid(i)      = find(dn>dnr(i),1,'first'); %first date after rain
end
timemat(timemat<0)=0;
dt     = dn(afid)-dnr'; %time since event, used in res.
LB     = dt; %lower bound on time - at least 1/2 time after first date.
UB     = maxT*ones(1,nr); %upper bound on time
tmod0  = startT*ones(1,nr); %starting time model

%open all files
if(isfield(fid,'perms'))
    fidi=rmfield(fidi,'perms'); %not using this here
end
[fidi,fido,online]=open_files(fidi,fido,nx,ny);


%% start
parpool(10)
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
     c0=nan(nc,nx);
     for i=1:nc
         [tmp,count1]=fread(fidi.c0(i).fid,nx,'real*4');
         if(count1>0)
             c0(i,1:count1)=tmp;
         end
     end
     c0=median(c0,1,'omitnan');
     c0(or(c0==-9999,or(isinf(c0),c0==0))) = NaN;
    
    goodid   = find(and(count>lowcountcutoff,c0>lowcorcutoff));
    
    %initialize arrays (necessary for parfor)
    mags    = nan(nr,length(goodid));times = mags; maglow = mags; maghig = mags; timerr = mags; fwd5 = mags; status=mags;
    shifts  = nan(1,length(goodid)); allres = shifts;
    
    parfor i=1:length(goodid)
        
        logd      = -log(dat(:,goodid(i)));
        deld      = logd(afid)-logd(befid);
        
        goodd     = isfinite(logd);
        weights   = diag(dsig(goodd,goodid(i)).^2); %data covariance matrix
        y         = logd(goodd);
        tmpmag    = zeros(nr,1);
        goodr     = deld>=-corcutoff*2; %throws away small changes
        ng        = sum(goodr);
        notdone   = ng>0; %start loop, as long as there are some "good" values
       
        
        %%%initialize values
        tmptime = NaN(nr,1);tmpstat=ones(nr,1);tstd = [];mstd = [];res = [];shift = 0; %mark tstd with 20 to see where ng cutoff occurrs
        tmpstat(~goodr)=2; %after inspection, these do seem to be all bad.
       
        while(notdone)
            x                   = timemat(goodd,goodr);
            [mod1,~,~,~,~,~,J0] = lsqnonlin('expfun',tmod0(goodr),LB(goodr),UB(goodr),options,x,y);

            jt                  = sum(abs(full(J0)));           
            [res,magst,shift,G] = expfun(mod1,x,y);
            
            tmpmag(goodr)         = magst;
            tmpmag(~goodr)        = 0;
            
            syntmp              = magst'.*exp(-dt(goodr)./mod1); %pred cor at dates after event
            
            badJ                = jt==0;
            toosmall            = syntmp<=corcutoff;       
            toobad              = or(badJ,toosmall);
            goodr(goodr)        = ~toobad;
            toobad              = sum(toobad);
            ng                  = sum(goodr);
            if(toobad==0)
                mod0                 = [mod1 magst' shift];
                J                    = jacobianest('expfun_all',mod0,x,y);
                J2                   = J'*J;
                J3                   = inv(J2);
                modcov               = J3*J'*weights*J*J3;
                modstd               = sqrt(diag(modcov));
                
                tstd                 = modstd(1:ng);
                mstd                 = modstd(ng+[1:ng]);
                
                badt                 = tstd>100;
                if(sum(badt))
                    notdone          = 1;
                    gid              = find(goodr);
                    goodr(gid(badt)) = false;
                    ng               = sum(goodr);
                else
                    notdone          = 0;
                    tmptime(goodr)   = mod1;
                    tmpmag(goodr)    = magst';
                    tmpstat(goodr)   = 3; %worked!
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

        tmpmag(tmpstat==2)=NaN; %in retrospect, these values are all bad.
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
        tmp(goodr)      = tmpmag(goodr)+magstd;
        maglow(:,i)     = tmp;
        tmp             = nan(nr,1);
        tmp(goodr)      = tmpmag(goodr)-magstd;
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
