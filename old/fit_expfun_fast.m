%parpool(5)
global divs allt
home=pwd;
dirs={'T54','T156','T47'};
dirs=dirs(1:2);
pol='_VV';

nt=length(dirs);
nd=0;
dates=[];
for i=1:length(dirs)
    tmp=dir([dirs{i} '/geo' pol '/rel*.cor.geo']);
    
    for j=1:length(tmp)
        t=tmp(j).name;
        dates(end+1).name=[tmp(j).folder '/' t];
        dates(end).dn=datenum(t(5:12),'yyyymmdd');
        dates(end).t=i;
    end
end
allt=[dates.t]';
nd=length(dates);
dn    = [dates.dn];

nx_geo    = 9361;
ny_geo    = 10297;


Gs=zeros(nd,nt);
for k=1:nd
    Gs(k,allt(k))=1;
end


rdir = ['results_TS' pol '/'];
rdate  = {'20150325','20150807','20160625','20170310','20170513','20170606'};
if(~exist(rdir,'dir'))
    mkdir(rdir)
end

dnr    = datenum(rdate,'yyyymmdd');
rdate  = rdate(dnr<max(dn));
dnr    = dnr(dnr<max(dn));
for i=1:length(dnr)
    rain(i).dnr=dnr(i);
    if(i<length(dnr))
        afid=find(and(dn>=dnr(i),dn<dnr(i+1)));
    else
        afid=find(dn>=dnr(i));
    end
    [jnk,sortid]=sort(dn(afid));
    rain(i).valid=afid(sortid);
end
divs=[min(dn)-1 dnr([1 6])'];

maxt    = 150;
mint    = 0;
LB      = mint*ones(length(dnr),1);
UB      = maxt*ones(length(dnr),1);
options = optimset('Display','off','TolFun',1e-3);

for i=1:nd
    fidr(i)=fopen(dates(i).name,'r');
end
for i=1:length(dirs)
    tmp=[dirs{i} '/geo' pol '/c0.cor.geo'];
    fid0(i)=fopen(tmp,'r');
end
online=0;
for i=1:length(dnr)
    fidmag(i)  = fopen([rdir rdate{i} '.mag0'],'w');
    fidmag10(i)  = fopen([rdir rdate{i} '.mag10'],'w');
    fidt(i)    = fopen([rdir rdate{i} '.time0'],'w');
end

fidres  = fopen([rdir 'resn0'],'w');


%start
for j=online+1:ny_geo
    j
    dat=zeros(nd,nx_geo);    
    for i=1:nd
        [tmp,count]=fread(fidr(i),nx_geo,'real*4');
        if(count>0)
            dat(i,1:count)=tmp;
        end
    end
    dat(dat==-9999)=NaN;
    dat(isinf(dat))=NaN;
    dat(dat==0)=NaN;
    count   = sum(isfinite(dat));
    goodid  = find(count>10);
   
    c0s=nan(length(dirs),nx_geo);
    for i=1:length(dirs)
        [tmp,count]=fread(fid0(i),nx_geo,'real*4');
        if(count>0)
            c0s(i,1:count)=tmp;
        end
    end
    c0s(c0s==-9999)=NaN;
    
    daterr  = 0.3083*(1-c0s(allt,:))-.2*(1-c0s(allt,:)).^2; %ndir x nx?
    
    doinv = false(length(dnr),length(goodid));
    mags  = zeros(length(dnr),length(goodid));
    times = zeros(length(dnr),length(goodid));
    resn  = nan(1,length(goodid));
    for i=1:length(dnr)
        valid  = [rain(i).valid];
        testd  = 1-dat(valid,goodid);
        teste  = daterr(valid,goodid);
        bigdat = testd>teste;%change this from 1.5
        goodn  = sum(bigdat,1);
        testd  = -log(dat(valid,goodid)); %logspace for rest of calc
        %x      = dn(valid)-rain(i).dnr; %time since rain event
        x=dn(valid)-min(dn(valid));
        maxx   = max(x);
        
        doinv(i,goodn>0)=true;
        
        %starting mag is just largest point
        mags(i,:) = max(testd,[],1,'omitnan');
        x         =  repmat(x',1,length(goodid));
        y         = testd./repmat(mags(i,:),length(valid),1);

        %only first point above noise = spike
        test=goodn==1;
        times(i,test)=6;  
        
        %some points above noise - fit exponential in logspace
        test=goodn>1;
        y(~bigdat)=NaN;
        slope = median(log(y)./x,1,'omitnan');
        slope = -1./slope;
        times(i,test)=slope(test);
        
        test=times(i,:)<0;
        doinv(i,test)=false;
        
        test=times(i,:)>maxx*2;
        times(i,test)=maxx*2;

        times(i,~doinv(i,:))=0;
        mags(i,~doinv(i,:))=0;
    end
                 
%     oldtimes=times;
%     oldmags=mags;
     for i=1:length(goodid)
        data     = dat(:,goodid(i));
        %de  = max(daterr(:,goodid(i)),0.3083*(1-data)-.2*(1-data).^2); %ndir x nx?
        
        g        = find(isfinite(data));
        data     = -log(data);
        tmpt     = times(:,i);
        %de=1./de;
        goodmag=doinv(:,i);
        if(sum(goodmag)>0)
            [newt,res]=lsqnonlin('expfunc_fast',tmpt(goodmag),LB(goodmag),UB(goodmag),options,dn(g),dnr(goodmag),data(g));
            [res2,synth2,test2]=expfunc_fast(newt,dn(g),dnr(goodmag),data(g));
            a=nan(length(dnr),1);
            b=a;
            a(goodmag)=(1-exp(-test2(1:sum(goodmag))));
            b(goodmag)=newt;
            mags(:,i)=a;
            times(:,i)=b;
            resn(i)=norm(res2);
        end
    end
    
    allresn             = nan(1,nx_geo);
    allmags             = nan(length(dnr),nx_geo);
    alltimes            = nan(length(dnr),nx_geo);
    allmags(:,goodid)   = mags;
    alltimes(:,goodid)  = times;
    allresn(goodid)     = resn;

  
    for i=1:length(dnr)
        fwrite(fidmag(i),allmags(i,:),'real*4');
        fwrite(fidt(i),alltimes(i,:),'real*4');
        fwrite(fidmag10(i),allmags(i,:).*exp(-10./alltimes(i,:)),'real*4');
    end
    fwrite(fidres,allresn,'real*4');
end
fclose('all');
