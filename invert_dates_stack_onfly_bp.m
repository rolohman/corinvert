parpool(8);
pol='_VV';
decide_ints_stack

%bin baseliens for later
load baselines.txt
bpr=baselines(id1)-baselines(id2);
abpr=abs(bpr);
% Set up fittype and options.
ft = fittype( 'a*x^2+b', 'independent', 'x', 'dependent', 'y' );

%nbin=20;
%[basebins,basebinedge]=discretize(abpr,nbin);
%basebinmid=basebinedge(1:end-1)+diff(basebinedge)/2;
bpw=200; %weights for baseline estimation

rx=rlooks;
ry=alooks;
rangevec=[1:newnx]*rx-floor(rx/2);
azvec=[1:newny]*ry-floor(ry/2);
alpha=50;

im=sqrt(-1);
windowtype=0;
switch windowtype
    case 0
        windx=[1/(rx+1):1/(rx+1):1 (1-1/(rx+1)):-1/(rx+1):1/(rx+1)];
        windy=[1/(ry+1):1/(ry+1):1 (1-1/(ry+1)):-1/(ry+1):1/(ry+1)];
    case 1
        windx=zeros(1,rx*2+1); windx(floor(rx/2):ceil(rx/2)+rx)=1;
        windy=zeros(1,ry*2+1); windy(floor(ry/2):ceil(ry/2)+ry)=1;
    case 2
        windx=exp(-(-rx:rx).^2/2/(rx/2)^2);
        windy=exp(-(-ry:ry).^2/2/(ry/2)^2);
end
windx=windx/sum(windx);
windy=windy/sum(windy);

ry=floor(length(windy)/2);


dni   = dn(1:nd-1)+diff(dn)/2;
%dni=dni(1:end-1);

intdt=diff(dn);

%open all files
fid=open_files(dates,'','invcor',pol);


adir     = ['results_dates_bp' pol '/'];
if(~exist(adir,'dir'))
    disp(['creating directory ' adir]);
    mkdir(adir)
end

if(exist([adir 'c0.cor'],'file'))
    fid0=fopen([adir 'c0.cor'],'r');
    for j=1:newny
        [tmp,count]=fread(fid0,newnx,'real*4');
        if(count<newnx)
            break
        end
    end
    fclose(fid0);
    
    online=j-1
    if(online<1)
        disp('c0.cor empty?')
        return
    end
    for i=1:nd
        fidr(i)=fopen([adir 'rel_' dates(i).name '.cor'],'a');
    end
    for i=1:nd-1
        fidp(i)=fopen([adir 'perm_' dates(i).name '_' dates(i+1).name '.cor'],'a');
    end
    fid0=fopen([adir 'c0.cor'],'a');
    fid1=fopen([adir 'resvar.cor'],'a');
    fidmin=fopen([adir 'cmin.cor'],'a');
    fidbp=fopen([adir 'bps'],'a');
    fid2=fopen([adir 'ovar.cor'],'a');
else
    online=0;
    for i=1:nd
        fidr(i)=fopen([adir 'rel_' dates(i).name '.cor'],'w');
    end
    for i=1:nd-1
        fidp(i)=fopen([adir 'perm_' dates(i).name '_' dates(i+1).name '.cor'],'w');
    end
    fid0=fopen([adir 'c0.cor'],'w');
    fid1=fopen([adir 'resvar.cor'],'w');
    fidmin=fopen([adir 'cmin.cor'],'w');
    fidbp=fopen([adir 'bps'],'w');
    fid2=fopen([adir 'ovar.cor'],'w');
end
cidl     = tril(ones(nd),-1)==1;
errslope = 0.3;
mncor    = errslope/(1+errslope);
mncorl   = log(mncor);

Gi0=zeros(ni,nd-1); %intervals, perm cor loss
for i=1:ni
    Gi0(i,id1(i):id2(i)-1)=1;
end
Gr0=zeros(ni,nd); %rel cor on dates
for i=1:ni
    Gr0(i,id1(i))=1;
    Gr0(i,id2(i))=-1;
end
Gr=Gr0(:,2:end);

for j=online+1:newny
    j
      if(j==1)
        readl=ry;
        startbit=0;
        starty=0;
    else
        readl=ry*2+1;
        starty=(j-1)*ry-ceil(ry/2);
        startbit=starty*nx*8; %end of line before buffer
    end
    
    slcs=nan(nx,ry*2+1,nd);
    for i=1:nd
        fseek(fid(i).in,startbit,-1);
        tmp=fread(fid(i).in,[nx*2,readl],'real*4');
        
        cpx=tmp(1:2:end,:)+im*tmp(2:2:end,:);
    
         
        if(j==newny)
            cpx(:,end+1:ry*2+1)=NaN;
        elseif(j==1)
            cpx=[nan(nx,ry+1) cpx];
        end
       
        slcs(:,:,i)=cpx;

    end
    count=0;
    cors=nan(ni,newnx);
    for i=1:nd-1    
        slc1=shiftdim(slcs(:,:,i));
        for k=i+1:nd
            count=count+1;
            slc2=shiftdim(slcs(:,:,k));
            a   = slc1.*conj(slc1);
            b   = slc2.*conj(slc2);
            c   = slc1.*conj(slc2);
            a   = windy*a';
            b   = windy*b';
            c   = windy*c';
            
            asum = conv(a,windx,'same');
            bsum = conv(b,windx,'same');
            csum = conv(c,windx,'same');
            cpx3 = csum./sqrt(asum.*bsum);
            cpx3 = cpx3(rangevec);
            
            sm   = abs(cpx3);
            sm(isnan(sm))=0;
            
            cors(count,:)=sm;
        end
    end
    
    
    good   = cors>0;
    cors(cors>1)=1;
    cors(~good)=NaN;
    ngood  = sum(good,1);
    goodid = find(ngood>=10);
    gcount = sum(cors>mncor,'omitnan');
    
    allm0  = NaN(nd*2+4,newnx); %c0,ct,cr,var,max,min,bp,dvar
    if(goodid)
        permbad=false(ni,newnx);
        doperm=find(and(gcount>0,gcount<ni*0.98));
        permbad(:,doperm)=countrows(id1,id2,cors(:,doperm),mncorl,mncor);
        
        t0=nan(1,length(goodid));
        t1=nan(nd-1,length(goodid));
        t2=nan(nd,length(goodid));
        t3=nan(1,length(goodid));
        t4=nan(1,length(goodid));
        t5=nan(1,length(goodid));
        t6=nan(1,length(goodid));
        tic
        parfor i=1:length(goodid)
          %  for i=1:length(goodid)
            
            data  = cors(:,goodid(i));
            d     = log(data);
            g     = isfinite(d);
            
            slopes=(d-max(d,[],'omitnan'))./abpr.^2;
            bigslopes=and(slopes<0,abpr>100);
            slopebound=min([-2e-8,max(slopes(bigslopes),[],'omitnan')]);
            dbound=min([mymax(d(g),100),max(d(g),[],'omitnan')-2e-8]);
            derrb = exp(bpw*d); %weights for bp fit.
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.StartPoint = [slopebound max(d(g))];
            opts.Lower=[slopebound dbound];
            opts.Upper=[0 0];
            opts.Weights = derrb(g);
            
            % Fit model to data.
            
            [fitresult, gof] = fit(abpr(g),d(g), ft, opts );
            if(and(fitresult.a<0,fitresult.b<0))
                bps=fitresult.a;
                d=d-bps*abpr.^2;
            else
                bps=0;
            end
            
            derr  = errslope*(1-data);
            dmin  = data-derr; dmin(dmin<0)=1e-10;
            dmax  = data+derr;
            
            dlmin = log(dmin); dlmax=log(dmax);
            derr  = dlmax-d;
               
            [gooddat,Grstat,Gistat]=find_bad(d-derr,diags,permbad(:,goodid(i)),Gr0,Gi0);
            d0    = d;
            c0    = mymax(d0(gooddat),100);
            
            cpmin = zeros(1,nd-1);
            w=exp(100*dlmin);
            dw=dlmin.*w;
            
            for k=1:nd-1
                id       = and(isfinite(dw),Gi0(:,k)==1);
                cpmin(k) = sum(dw(id))/sum(w(id));
            end
            cpmin(~Gistat) = d(diags(~Gistat))-c0;
            cpmin          = min(0,cpmin);
            
            [cp1]          = est_ct(d,Gi0,cpmin);
            cp1(~Gistat)   = min(0,d(diags(~Gistat))-c0);
            
            synp           = Gi0*cp1';
            synp(~gooddat) = NaN;
            c0             = min(0,mymax(d0-synp,25));
            cres           = d-synp-c0;
            cr1            = est_cr(cres,exp(synp),nd,cidl);
            cr1(~Grstat)   = NaN;
            [cshifts1,cp2] = flatten_front(cr1,25,cp1,cpmin);
            cp2(~Gistat)   = min(0,d(diags(~Gistat))-c0);
            
            synp           = Gi0*cp2';
            synp(~gooddat) = NaN;
            cres           = d-synp-c0;
            cr2            = est_cr(cres,exp(synp),nd,cidl);
            cr2(~Grstat)   = NaN;
            [cshifts2,cp]  = flatten_back(cr2,25,cp2,cpmin);
            cp(~Gistat)    = min(0,d(diags(~Gistat))-c0);
            
            synp           = Gi0*cp';
            synp(~gooddat) = NaN;
            cres           = d-synp-c0;
            cr             = est_cr(cres,exp(synp),nd,cidl);
            %final, no more inverts
            [cshifts,cp] = flatten_front(cr,25,cp,cpmin);
            cr=cr-cshifts;
            cp(~Gistat)   = min(0,d(diags(~Gistat))-c0);
            [cshifts,cp]  = flatten_back(cr,25,cp,cpmin);
            cr=cr-cshifts;
            cp(~Gistat)    = min(0,d(diags(~Gistat))-c0);
            
            rf             = isfinite(cr);
            res            = d-Gi0*cp'+abs(Gr0(:,rf)*cr(rf)')-c0;           
            cr             = cr-mymax(cr,50);
           
            t0(i)   = exp(c0);
            t1(:,i) = exp(cp');
            t2(:,i) = exp(cr');
            t3(i)   = var(res(gooddat),[],'omitnan');
            t4(i)   = 1-mymax(1-data,100);
            t5(i)   = bps;
            t6(i)   = var(d0(gooddat));
        end
        tb=toc;
        disp((newny-j)*tb/60/60)
        %bookkeeping related to parallelization
        
        allm0(1,goodid)         = t0; %c0
        allm0(2:nd,goodid)      = t1; %cp
        allm0(nd+1:nd*2,goodid) = t2; %cr
        allm0(end-3,goodid)     = t3; %var
        allm0(end-2,goodid)     = t4; %min
        allm0(end-1,goodid)     = t5; %bps
        allm0(end,goodid)       = t6; %origvar
    end
    
    %write output files for this line
    fwrite(fid0,allm0(1,:),'real*4');
    for i=1:nd-1
        fwrite(fidp(i),allm0(1+i,:),'real*4');
    end
    for i=1:nd
        fwrite(fidr(i),allm0(nd+i,:),'real*4');
    end
    fwrite(fid1,allm0(end-3,:),'real*4');
    fwrite(fidmin,allm0(end-2,:),'real*4');
    fwrite(fidbp,allm0(end-1,:),'real*4');
    fwrite(fid2,allm0(end,:),'real*4');
end

fclose('all');


