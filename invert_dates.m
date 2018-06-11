parpool(6);
decide_ints
%now assume all ints made
ints=[];
for j=1:nd-1
    for i=j+1:nd
        ints(end+1).id1=j;
        ints(end).id2=i;
    end
end
id1    = [ints.id1]';
id2    = [ints.id2]';
diags  = find(id2==id1+1);
ni     = length(ints);
nd     = length(dates);
bp     = [dates.bp];

dn     = [dates.dn];
dn1    = dn(id1);
dn2    = dn(id2);
intdt  = diff(dn); %intervals between dates (not interferograms)
adir   = 'results_dates/';
rlooks = 7;
alooks = 3;
n      = rlooks*alooks;
newnx  = floor(nx/rlooks)
newny  = floor(ny/alooks);
cidl   = tril(ones(nd),-1)==1;
errslope = 0.3;
mncor  = errslope/(1+errslope);
mncorl = log(mncor);


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
%Gi=Gi0(:,1:end-1);
%Gp=[ones(ni,1) Gi];

for i=1:ni
    cordir=(['cordir/' dates(id1(i)).name '/']);
    intcordir=[cordir dates(id2(i)).name '/'];
    corfile_small=[intcordir dates(id1(i)).name '_' dates(id2(i)).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
    if(exist(corfile_small,'file'))
        fidi(i)=fopen(corfile_small,'r');
    else
        disp(['where is ' corfile_small '?'])
    end
end


if(~exist(adir,'dir'))
    mkdir(adir)
end

%find last line written
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
    fid0=fopen([adir 'c0.cor'],'a');
    fid1=fopen([adir 'rms.cor'],'a');
    fidmin=fopen([adir 'cmin.cor'],'w');
    
    for i=1:nd
        fidr(i)=fopen([adir 'rel_' dates(i).name '.cor'],'a');
    end
    for i=1:nd-1
        fidp(i)=fopen([adir 'perm_' dates(i).name '_' dates(i+1).name '.cor'],'a');
    end
    for i=1:ni
        fseek(fidi(i),online*newnx*4,-1);
    end
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


for j=online+1:newny
    j
    cors=zeros(ni,newnx);
    
    for i=1:ni
        [tmp,count]=fread(fidi(i),newnx,'real*4');
        if(count>0)
            cors(i,1:count)=tmp;
        end
    end
    
    good   = cors>0;
    cors(cors>1)=1;
    cors(~good)=NaN;
    ngood  = sum(good,1);
    goodid = find(ngood>=100);
    gcount = sum(cors>mncor,'omitnan');
    
    allm0  = NaN(nd*2+2,newnx); %c0,ct,cr,rms,max,min
    if(goodid)
        permbad=false(ni,newnx);
        doperm=find(and(gcount>0,gcount<ni*0.98));
        
        permbad(:,doperm)=countrows(id1,id2,cors(:,doperm),mncorl,mncor);
        t0=nan(1,length(goodid)); %c0 cmax cmin
        t1=nan(nd-1,length(goodid));
        t2=nan(nd,length(goodid));
        t3=nan(1,length(goodid));
        t4=nan(1,length(goodid));
        alpha=100;

       tic;
        parfor i=1:length(goodid)
          
            data  = cors(:,goodid(i));
            
            derr  = errslope*(1-data);
            dmin=data-derr; dmin(dmin<0)=1e-10;
            dmax=data+derr;
            
            d=log(data); dlmin=log(dmin); dlmax=log(dmax);
            derr=dlmax-d;
            
            
            [gooddat,Grstat,Gistat]=find_bad(d-derr,diags,permbad(:,goodid(i)),Gr0,Gi0);
            
            %c0=max(dlmin(gooddat));
            c0=mymax(d(gooddat),100);
            cpmin=zeros(1,nd-1);
            for k=1:nd-1
                id=Gi0(:,k)==1;
                cpmin(k)=mymax(dlmin(id),100);
            end
            cpmin(~Gistat)=d(diags(~Gistat));
            cpmin=min(0,cpmin-c0);
            [cp,c0]=est_ct_c0_new(d,Gi0,cpmin);
            cp(~Gistat)=d(diags(~Gistat));
            
            synp=Gi0*cp';
            cres=d-c0-synp;
            
            cr=est_cr(cres,exp(synp),nd,cidl);
            [cshifts,cp2]=flatten(cr,alpha,cp,cpmin);
            cr=cr-cshifts;
            cp=cp2;
            %[cp,cr]=est_all_new(c0,d,gooddat,Gi0,Gr0,Gistat,Grstat,diags,id1,id2,nd,cpmin);
            res=d-Gi0*cp'+abs(Gr0*cr')-c0;
            
            cr=cr-mymax(cr,25);
            
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


