parpool(4);
decide_ints

adir     = 'results_dates/';
open_files

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


