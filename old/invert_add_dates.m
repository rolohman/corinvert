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
        windx=ones(1,rx*2+1);
        windy=ones(1,ry*2+1);
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

%open all slcs
clear dates fidi
for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
    fidi(i)          = fopen(dates(i).slc,'r');
end

adir     = ['results_dates' pol '/'];
if(~exist(adir,'dir'))
    disp([adir ' does not exist']);
end

dirc0=dir([adir 'c0.cor']);
bytlen=dirc0.bytes;

for i=1:nd
    outfiles(i).rel=[adir 'rel_' dates(i).name '.cor'];
    if(exist(outfiles(i).rel,'file'))
        dirtemp=dir(outfiles(i).rel);
        if(dirtemp.bytes<bytlen)
            outfiles(i).reldone=1;
            outfiles(i).online=dirtemp.bytes/newnx/4;
        else
            outfiles(i).reldone=2;
        end
    else
        outfiles(i).reldone=0;
    end
end
done0=find([outfiles.reldone]==0);
done1=find([outfiles.reldone]==1);
done2=find([outfiles.reldone]==2);
if(and(length(done0)>0,length(done1)>0))
    disp('shouldnt have empty files and nonexistent files')
    return
elseif(and(length(done0)==0,length(done1)==0))
    disp('nothing to do')
end

fid0=fopen([adir 'c0.cor'],'r');

if(length(done1)>0)
    online=outfiles(done1(1)).online;
    disp(['restarting inversion from line ' num2str(online)])
    done1=[done1(1)-1 done1];
    dates=dates(done1);
    fidi=fidi(done1);
    dn=dn(done1);
    nd=length(dn);
    
    fseek(fid0,online*newnx*4,-1)
    for i=1:nd
        fidr(i)=fopen([adir 'rel_' dates(i).name '.cor'],'a');
    end
    for i=1:nd-1
        fidp(i)=fopen([adir 'perm_' dates(i).name '_' dates(i+1).name '.cor'],'a');
    end

else
    online=0;
    done0=[done2(end-1:end) done0];
    olddates=dates(done2(1:end-2));
    fidold=fidi(done2(1:end-2));
    
    dates=dates(done0);
    fidi=fidi(done0);

    dntot=dn;
    dn=dn(done0);
    nd=length(dn);
    
    movefile([adir 'rel_' dates(1).name '.cor'],[adir 'rel_' dates(1).name '_old.cor']);
    movefile([adir 'rel_' dates(2).name '.cor'],[adir 'rel_' dates(2).name '_old.cor']);
    movefile([adir 'perm_' dates(1).name '_' dates(2).name '.cor'],[adir 'perm_' dates(1).name '_' dates(2).name '_old.cor']);
    
    for i=1:length(olddates)
        fidrold(i)=fopen([adir 'rel_' olddates(i).name '.cor'],'r');
    end
    for i=1:nd
        fidr(i)=fopen([adir 'rel_' dates(i).name '.cor'],'w');
    end
    for i=1:nd-1
        fidp(i)=fopen([adir 'perm_' dates(i).name '_' dates(i+1).name '.cor'],'w');
    end
end

ints=[];
for j=1:nd-1
    for i=j+1:nd %make all
        ints(end+1).id1=j;
        ints(end).id2=i;
    end
end
id1    = [ints.id1]';
id2    = [ints.id2]';
diags  = find(id2==id1+1);
ni     = length(ints);

dn1    = dn(id1);
dn2    = dn(id2);

dni   = dn(1:nd-1)+diff(dn)/2;
%dni=dni(1:end-1);

intdt=diff(dn);



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
    allc0=fread(fid0,newnx,'real*4');
    for i=1:length(olddates)
       alloldcr(i,:)=fread(fidrold(i),newnx,'real*4');
    end
  
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
    
    allm0  = NaN(nd*2,newnx); %c0,ct,cr,rms,max,min
    if(goodid)
        permbad=false(ni,newnx);
        doperm=find(gcount>0);
        permbad(:,doperm)=countrows(id1,id2,cors(:,doperm),mncorl,mncor);
        
        t1=nan(nd-1,length(goodid));
        t2=nan(nd,length(goodid));
    
        tic
        for i=1:length(goodid)
            
            data  = cors(:,goodid(i));
            
            derr  = errslope*(1-data);
            dmin=data-derr; dmin(dmin<0)=1e-10;
            dmax=data+derr;
            
            d=log(data); dlmin=log(dmin); dlmax=log(dmax);
            derr=dlmax-d;
               
            [gooddat,Grstat,Gistat]=find_bad(d-derr,diags,permbad(:,goodid(i)),Gr0,Gi0);
            
            c0    = log(allc0(goodid(i)));
            oldcr = log(alloldcr(:,goodid(i)));
            zcr=0*oldcr;
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
            cres           = d-synp-c0;
            cr1            = est_cr(cres,exp(synp),nd,cidl);
            cr1(~Grstat)   = NaN;
            [cshifts1,cp2] = flatten_front([oldcr' cr1],25,[zcr' cp1],[zcr' cpmin]);
            cp2=cp2(end-nd+2:end);
            cp2(~Gistat)   = min(0,d(diags(~Gistat))-c0);
            
            synp           = Gi0*cp2';
            synp(~gooddat) = NaN;
            cres           = d-synp-c0;
            cr2            = est_cr(cres,exp(synp),nd,cidl);
            cr2(~Grstat)   = NaN;
            [cshifts2,cp]  = flatten_back([oldcr' cr2],25,[zcr' cp2],[zcr' cpmin]);
            cp=cp(end-nd+2:end);
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
            
            t1(:,i) = exp(cp');
            t2(:,i) = exp(cr');
 
            
        end
        tb=toc;
        disp((newny-j)*tb/60/60)
        %bookkeeping related to parallelization
        
        allm0(2:nd,goodid)      = t1; %cp
        allm0(nd+1:nd*2,goodid) = t2; %cr
       end
    
    %write output files for this line
    for i=1:nd-1
        fwrite(fidp(i),allm0(1+i,:),'real*4');
    end
    for i=1:nd
        fwrite(fidr(i),allm0(nd+i,:),'real*4');
    end
 
end

fclose('all');


