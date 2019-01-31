function [allcors,allphs]= plot_int_onfly(xpt,ypt,plotflag)

decide_ints_stack
slcdir=['merged/SLC' pol '/'];
rx=rlooks;
ry=alooks;

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


%matrix that maps ints to a grid
cids=[];
for i=1:nd
    cids=[cids (nd)*(i-1)+[i:nd-1]];
end

%open all slcs
for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
    fidi(i)          = fopen(dates(i).slc,'r');
end

for k=1:length(xpt)
    %start/stop points   
    startx=max(1,ceil(xpt(k)*rx-(rx/2)-rx));
    stopx=min(nx,ceil(xpt(k)*rx-(rx/2)+rx));
    if(startx==1)
        bx=ceil(rx/2)+1:rx*2+1;
    elseif(stopx==nx)
        bx=1:rx*2+1-ceil(rx/2);
    else
        bx=1:rx*2+1;
    end
    if(ypt(k)==1)
        readl=ry;
        startbit=0;
        starty=0;
    else
        readl=ry*2+1;
        starty=(ypt(k)-1)*ry-ceil(ry/2);
        startbit=starty*nx*8; %end of line before buffer
    end
    
    %calculate start/stop points
    slcs=nan(nd,rx*2+1,ry*2+1);
    for i=1:nd
        fseek(fidi(i),startbit,-1);
        tmp=fread(fidi(i),[nx*2,readl],'real*4');
        
        cpx=tmp(1:2:end,:)+im*tmp(2:2:end,:);
        cpx=cpx(startx:stopx,:);
         
        if(ypt(k)==newny)
            cpx(:,end+1:ry*2+1)=NaN;
        elseif(ypt(k)==1)
            cpx=[nan(stopx-startx+1,ry+1) cpx];
        end
       
        slcs(i,bx,:)=cpx;

    end
    for i=1:ni
        
        cpx1=squeeze(slcs(id1(i),:,:));
        cpx2=squeeze(slcs(id2(i),:,:));
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
        
phs(i)=angle(cpx3(rx+1));
        cors(i)=sm(rx+1);
    end
    
    
    good        = cors>0;
    cors(~good) = NaN;
    ngood       = sum(good,1);
    allcors(:,k)=cors;
allphs(:,k)=phs;
    if(plotflag)
    
    badi=[];
    jnk=nan(nd);
if(plotflag)
figure('name',[num2str(xpt(k)) ' ' num2str(ypt(k))])
    colormap('jet')
    jnk(cids)=cors;
    pcolor(dn,dn,jnk'),shading flat,set(gca,'ydir','reverse');
    hold on
    fixplot
    
    title('data')
    end
end
end
fclose('all');
