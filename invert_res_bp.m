%parpool(8);
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

for i=1:nd
    fidr(i)=fopen([adir 'rel_' dates(i).name '.cor'],'r');
end
for i=1:nd-1
    fidp(i)=fopen([adir 'perm_' dates(i).name '_' dates(i+1).name '.cor'],'r');
end
fid0=fopen([adir 'c0.cor'],'r');
fidb=fopen([adir 'bpslope.cor'],'w');
fidbe=fopen([adir 'bpslopeerr.cor'],'w');
fidbo=fopen([adir 'bpoff.cor'],'w');

cidl     = tril(ones(nd),-1)==1;
bp=load('baselines.txt');
bpr=bp(id1)-bp(id2);

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

for j=1:newny
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
        fseek(fidi(i),startbit,-1);
        tmp=fread(fidi(i),[nx*2,readl],'real*4');
        
        cpx=tmp(1:2:end,:)+im*tmp(2:2:end,:);
    
         
        if(j==newny)
            cpx(:,end+1:ry*2+1)=NaN;
        elseif(j==1)
            cpx=[nan(nx,ry+1) cpx];
        end
       
        slcs(:,:,i)=cpx;

    end
    
     cr=zeros(nd,newnx);
    cp=zeros(nd-1,newnx);
    
    for i=1:nd
        cr(i,:)=fread(fidr(i),newnx,'real*4');
    end
    for i=1:nd-1
        cp(i,:)=fread(fidp(i),newnx,'real*4');
    end
    c0=fread(fid0,newnx,'real*4');
 
  
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
    cr(isnan(cr))=1;
    c0(isnan(c0))=1;
    cp(isnan(cp))=1;
    
 
    cors(cors==0)=NaN;
    for i=1:newnx
        d=log(cors(:,i));
        good=sum(isfinite(d));
        if(good>50)
            
            syn=Gi0*log(cp(:,i))-abs(Gr0*log(cr(:,i)))+log(c0(i));
            
            res            = d-syn;
            [b,stats]=robustfit(abs(bpr),res);
            fwrite(fidb,b(2),'real*4');
            fwrite(fidbo,b(1),'real*4');
            fwrite(fidbe,stats.se(2),'real*4');
        else
            fwrite(fidb,nan(1),'real*4');
            fwrite(fidbo,nan(1),'real*4');
            fwrite(fidbe,nan(1),'real*4');
            
        end
        
    end
   
end

fclose('all');


