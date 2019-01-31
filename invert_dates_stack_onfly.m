parpool(8);
pol='_VV';
decide_ints_stack

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

z   = zeros(1,nx);
rea = zeros(ry*2+1,nx);
ima = zeros(ry*2+1,nx);

dni   = dn(1:nd-1)+diff(dn)/2;
%dni=dni(1:end-1);

intdt=diff(dn);


%open all slcs
for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
    fidi(i)          = fopen(dates(i).slc,'r');
end

adir     = ['results_dates' pol '/'];
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
    fid1=fopen([adir 'rms.cor'],'a');
    fidmin=fopen([adir 'cmin.cor'],'a');
else
    online=0;
    for i=1:nd
        fidr(i)=fopen([adir 'rel_' dates(i).name '.cor'],'w');
    end
    for i=1:nd-1
        fidp(i)=fopen([adir 'perm_' dates(i).name '_' dates(i+1).name '.cor'],'w');
    end
    fid0=fopen([adir 'c0.cor'],'w');
    fid1=fopen([adir 'rms.cor'],'w');
    fidmin=fopen([adir 'cmin.cor'],'w');
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
    
    slcs=nan(nd,nx,ry*2+1);
    for i=1:nd
        fseek(fidi(i),startbit,-1);
        tmp=fread(fidi(i),[nx*2,readl],'real*4');
        
        cpx=tmp(1:2:end,:)+im*tmp(2:2:end,:);
    
         
        if(j==newny)
            cpx(:,end+1:ry*2+1)=NaN;
        elseif(j==1)
            cpx=[nan(nx,ry+1) cpx];
        end
       
        slcs(i,:,:)=cpx;

    end
  
    for i=1:ni
        
        cpx1=shiftdim(slcs(id1(i),:,:));
        cpx2=shiftdim(slcs(id2(i),:,:));
        a1=abs(cpx1);
        a2=abs(cpx2);
        cpx=cpx1.*conj(cpx2);
        am=sqrt(a1.*a2);
        goodid=am>0;
        cpx(goodid)=cpx(goodid)./am(goodid);
        
        rea=real(cpx)';
        ima=imag(cpx)';
        
        mag  = sqrt(rea.^2+ima.^2);
        
        r2   = rea;
        i2   = ima;
        r2   = windy*r2;
        i2   = windy*i2;
        m2   = windy*mag;
        
        r2(isnan(r2))=0;
        i2(isnan(i2))=0;
        m2(isnan(m2))=1;
        rsum = conv(r2,windx,'same');
        isum = conv(i2,windx,'same');
        msum = conv(m2,windx,'same');
        cpx3 = rsum+im*isum;
        cpx3 = cpx3./msum;
        sm   = abs(cpx3);
        sm(isnan(sm))=0;
        
        cors(i,:)=sm(rangevec);
    end
   
    
    good   = cors>0;
    cors(cors>1)=1;
    cors(~good)=NaN;
    ngood  = sum(good,1);
    goodid = find(ngood>=10);
    gcount = sum(cors>mncor,'omitnan');
    
    allm0  = NaN(nd*2+2,newnx); %c0,ct,cr,rms,max,min
    if(goodid)
        permbad=false(ni,newnx);
        doperm=find(and(gcount>0,gcount<ni*0.98));
        permbad(:,doperm)=countrows(id1,id2,cors(:,doperm),mncorl,mncor);
        
        t0=nan(1,length(goodid)); 
        t1=nan(nd-1,length(goodid));
        t2=nan(nd,length(goodid));
        t3=nan(1,length(goodid));
        t4=nan(1,length(goodid));
    
        tic
        parfor i=1:length(goodid)
            
            data  = cors(:,goodid(i));
            
            derr  = errslope*(1-data);
            dmin=data-derr; dmin(dmin<0)=1e-10;
            dmax=data+derr;
            
            d=log(data); dlmin=log(dmin); dlmax=log(dmax);
            derr=dlmax-d;
               
            [gooddat,Grstat,Gistat]=find_bad(d-derr,diags,permbad(:,goodid(i)),Gr0,Gi0);
            
            c0    = mymax(d(gooddat),100);
            
            cpmin = zeros(1,nd-1);
            for k=1:nd-1
                id       = Gi0(:,k)==1;
                cpmin(k) = mymax(dlmin(id),100);
            end
            cpmin(~Gistat) = d(diags(~Gistat))-c0;
            cpmin          = min(0,cpmin);
            [cp1]          = est_ct(d,Gi0,cpmin);
            cp1(~Gistat)   = min(0,d(diags(~Gistat))-c0);
            
            synp           = Gi0*cp1';
            synp(~gooddat) = NaN;
            c0             = min(0,mymax(d-synp,25));
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
            
            res            = d-Gi0*cp'+abs(Gr0*cr')-c0;           
            cr             = cr-mymax(cr,50);
            
            t0(i)   = exp(c0);
            t1(:,i) = exp(cp');
            t2(:,i) = exp(cr');
            t3(i)   = sqrt(mean((res(gooddat)).^2));
            t4(i)   = 1-mymax(1-data,100);
            
        end
        tb=toc;
        disp((newny-j)*tb/60/60)
        %bookkeeping related to parallelization
        
        allm0(1,goodid)         = t0; %c0
        allm0(2:nd,goodid)      = t1; %cp
        allm0(nd+1:nd*2,goodid) = t2; %cr
        allm0(end-1,goodid)     = t3; %rms
        allm0(end,goodid)       = t4; %min
    end
    
    %write output files for this line
    fwrite(fid0,allm0(1,:),'real*4');
    for i=1:nd-1
        fwrite(fidp(i),allm0(1+i,:),'real*4');
    end
    for i=1:nd
        fwrite(fidr(i),allm0(nd+i,:),'real*4');
    end
    fwrite(fid1,allm0(end-1,:),'real*4');
    fwrite(fidmin,allm0(end,:),'real*4');

end

fclose('all');


