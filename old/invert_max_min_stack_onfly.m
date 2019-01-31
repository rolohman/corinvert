pol='_VV';
decide_ints_stack

rx=rlooks;
ry=alooks;
rangevec=[1:newnx]*rx-floor(rx/2);
azvec=[1:newny]*ry-floor(ry/2);

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

dni   = dn(1:nd-1)+diff(dn)/2;
%dni=dni(1:end-1);

intdt=diff(dn);


%open all slcs
for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
    fidi(i)          = fopen(dates(i).slc,'r');
end

adir     = 'results_dates/';
if(~exist(adir,'dir'))
    disp(['creating directory ' adir]);
    mkdir(adir)
end

alpha=50;

fid1=fopen([adir 'cmax_robust.cor'],'w');
fid2=fopen([adir 'cmax.cor'],'w');
fid3=fopen([adir 'cmin_robust.cor'],'w');
fid4=fopen([adir 'cmin.cor'],'w');

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
   
    allm  = NaN(4,newnx); %cmax_r cmax cmin_r cmin
    if(goodid)
        allm(1,:)=mymax_mat(cors,alpha);
        allm(2,:)=max(cors,[],1,'omitnan');
        allm(3,:)=1-mymax_mat(1-cors,alpha);
        allm(4,:)=min(cors,[],1,'omitnan');
    end
    
    %write output files for this line
    fwrite(fid1,allm(1,:),'real*4');
    fwrite(fid2,allm(2,:),'real*4');
    fwrite(fid3,allm(3,:),'real*4');
    fwrite(fid4,allm(4,:),'real*4');
end

fclose('all');


