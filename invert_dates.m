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
ni     = length(ints);
nd     = length(dates);
bp     = [dates.bp];

dn     = [dates.dn];
dn1    = dn(id1);
dn2    = dn(id2);
intdt  = diff(dn); %intervals between dates (not interferograms)
adir   = 'results_dates_newstuff/';
rlooks = 7;
alooks = 3;
n      = rlooks*alooks;
newnx  = floor(nx/rlooks)
newny  = floor(ny/alooks);
cidl   = tril(ones(nd),-1)==1;
options=optimset('Display','Off','TolFun',1e-5);
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
Gi=Gi0;
Gp=[ones(ni,1) Gi];

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
    
    allm0  = NaN(nd*2+3,newnx); %c0,ct,cr,rms,max,min
    if(goodid)
        t0=nan(1,length(goodid)); %c0 cmax cmin
        t1=nan(nd-1,length(goodid));
        t2=nan(nd,length(goodid));
        t3=nan(1,length(goodid));
        t4=nan(1,length(goodid));
        
       parfor i=1:length(goodid)
            data  = cors(:,goodid(i));
            goodi = isfinite(data);
            ng    = sum(goodi);
         
            
            [ct_est0,c0_est0,synth,cres]=est_ct_c0_smoothmax(data,id1,id2,Gi);
           
            c0_est0 = log(c0_est0);
            ct_est0 = log(ct_est0');
      
            cr0   = est_cr(cres,synth,nd,cidl,0);
            cr0   = log(cr0);
 
            tres   = log(data)-(-abs(Gr0*cr0')); 
            mod    = -lsqnonneg(-Gp(goodi,:),tres(goodi));
            
            ct_est1 = mod(2:end);
            
            sp     = Gi*ct_est1;
            cres   = log(data)-c0_est0-sp;
            cr1    = est_cr(exp(cres),exp(sp),nd,cidl,0);
            cr1    = log(cr1);
            
            finalsynth       = c0_est0+Gi*ct_est1-abs(Gr0*cr1');
            
            finalmod         = cr1;
            finalmod         = exp(finalmod);
                
            t0(i)   = exp(c0_est0);
            t1(:,i) = exp(ct_est1);
            t2(:,i) = finalmod; 
            t3(i)   = sqrt(mean((data(goodi)-finalsynth(goodi)).^2));
            t4(i)   = 1-mymax(1-data,100);
       end
      
        %bookkeeping related to parallelization
        allm0(1,goodid)         = t0;
        allm0(2:nd,goodid)      = t1;
        allm0(nd+1:nd*2,goodid) = t2;
        allm0(end-1,goodid)     = t3;
        allm0(end,goodid)       = t4;
    end
    
    %write output files for this line
    fwrite(fid0,allm0(1,:),'real*4');
    fwrite(fid1,allm0(end-1,:),'real*4');
    fwrite(fidmin,allm0(end,:),'real*4');
    for i=1:nd-1
        fwrite(fidp(i),allm0(1+i,:),'real*4');
    end
    for i=1:nd
        fwrite(fidr(i),allm0(nd+i,:),'real*4');
    end
end

fclose('all');


