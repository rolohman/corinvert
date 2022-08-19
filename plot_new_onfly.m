function [output]= plot_new_onfly(xpt,ypt,plotflag,slcflag)
%out=plot_all_onflyb(2725*rx+348-1,2900*ry+101-1+169); gives you the pixel
%from extract_chunk_2725_2775_2900_3000 at 348,101 after removal of 1:169
%columns

if(~exist('slcflag','var'))
    slcflag=0;
end
if(~exist('plotflag','var'))
    plotflag=0;
end

lonfile='merged/geom_master/lon.rdr.4alks_15rlks.full';
latfile='merged/geom_master/lat.rdr.4alks_15rlks.full';

%window to pull out around pt.
nrx=5;
nry=5;

pol='_VV';

clear dates
decide_ints_stack

output.dn=dn;
output.id1=id1;
output.id2=id2;

rx=rlooks;
ry=alooks;

windowtype=2;
switch windowtype
    case 0
        windx=[1/(nrx+1):1/(nrx+1):1 (1-1/(nrx+1)):-1/(nrx+1):1/(nrx+1)];
        windy=[1/(nry+1):1/(nry+1):1 (1-1/(nry+1)):-1/(nry+1):1/(nry+1)];
    case 1
        windx=zeros(1,nrx*2+1); windx(floor(nrx/2):ceil(nrx/2)+nrx)=1;
        windy=zeros(1,nry*2+1); windy(floor(nry/2):ceil(nry/2)+nry)=1;
    case 2
        windx=exp(-(-nrx:nrx).^2/2/(nrx/2)^2);
        windy=exp(-(-nry:nry).^2/2/(nry/2)^2);
end
windx=windx/sum(windx);
windy=windy/sum(windy);
wind = windx'*windy;
windn(1,:,:) = wind;
wgts         = conv2(ones(nry*2+1,nrx*2+1),wind,'same');
wgtsn(1,:,:) = wgts; %for 3d arrays.
slc1 = zeros(nry*2+1,nx);
slc2 = slc1;

intid     = nchoosek(1:nd,2);

%open all slcs
for i=1:nd
    dates(i).name    = files(i).name(1:8);
    dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
    fidi(i)          = fopen(dates(i).slc,'r');
end
%check if inversion is done (assume c0 is typical)
tdir       = ['results_dates' pol '/'];
c0filename = [tdir 'c0.cor'];
d          = dir(c0filename);

if(exist(lonfile,'file'))
    fid=fopen(lonfile,'r');
    fseek(fid,(newnx*(floor(ypt/ry)-1)+floor(xpt/rx)-1)*8,-1);
    output.lon=fread(fid,1,'real*8');
    
    fid=fopen(latfile,'r');
    fseek(fid,(newnx*(floor(ypt/ry)-1)+floor(xpt/rx)-1)*8,-1);
    output.lat=fread(fid,1,'real*8');
end
output.x=xpt;
output.y=ypt;

rels=zeros(nd,1);
modp=zeros(nd-1,1);

fid0=fopen([tdir '/c0.cor'],'r');
fseek(fid0,(newnx*(floor(ypt/ry)-1)+floor(xpt/rx)-1)*4,-1);
c0=fread(fid0,1,'real*4');
fclose(fid0);

for i=1:nd
    fidr(i)=fopen([tdir '/rel_' dates(i).name '.cor'],'r');
end
for i=1:nd
    if(fidr(i)>0)
        fseek(fidr(i),(newnx*(floor(ypt/ry)-1)+floor(xpt/rx)-1)*4,-1);
        rels(i)=fread(fidr(i),1,'real*4');
        fclose(fidr(i));
    else
        rels(i)=1;
    end
end
rels(isnan(rels))=1;

for i=1:nd-1
    fidp(i)=fopen([tdir '/perm_' dates(i).name '_' dates(i+1).name '.cor'],'r');
end
for i=1:nd-1
    if(fidp(i)>0)
        fseek(fidp(i),(newnx*(floor(ypt/ry)-1)+floor(xpt/rx)-1)*4,-1);
        modp(i)=fread(fidp(i),1,'real*4');
        fclose(fidp(i));
    else
        modp(i)=1;
    end
    
end
modp(isnan(modp))=1;
    
output.c0   = c0;
output.rels = rels;
output.modp = modp;

%now get slcs
clear cors phs hp
%start/stop points
startx=max(1,ceil(xpt-nrx));
stopx=min(nx,ceil(xpt+nrx));
starty=max(1,ceil(ypt-nry));
stopy =min(ny,ceil(ypt+nry));
startbit=(starty-1)*nx*8; %end of line before buffer

%calculate start/stop points
slcs=nan(nd,stopx-startx+1,stopy-starty+1);
for i=1:nd
    fseek(fidi(i),startbit,-1);
    tmp=fread(fidi(i),[nx*2,stopy-starty+1],'real*4');
    
    cpx=tmp(1:2:end,:)+1j*tmp(2:2:end,:);
    slcs(i,:,:)=cpx(startx:stopx,:);
end
fclose('all');

mags=reshape(slcs,nd,size(slcs,2)*size(slcs,3));
mags=abs(mags);
mags=median(mags,2,'omitnan');
output.mags=mags;
slcs=slcs./abs(slcs);
if(slcflag)
    output.slcs=slcs;
end

%calc cor, hp
allints   = slcs(intid(:,1),:,:).*conj(slcs(intid(:,2),:,:));
allintf   = convn(allints,windn,'same')./wgtsn;
    
allcor    = abs(allintf);
allintf   = allintf./allcor;
allinthp  = allints.*conj(allintf);
allinthp  = allinthp./abs(allinthp);
triplets0 = trip_nm(squeeze(allintf(:,nrx+1,nry+1)),intid,1);


output.trip0=triplets0;

output.hp   = allinthp(:,nrx+1,nry+1);
output.lp   = allintf(:,nrx+1,nry+1);
%output.phs  = allints(:,nrx+1,nry+1);
output.cors = allcor(:,nrx+1,nry+1);

if(c0>=0.5)
    disp('good cor')
    output.raincount = sum(rels<0.8,'omitnan');
    disp([num2str(output.raincount) ' rainy dates'])
    if(output.raincount>0)
        disp('doing correction')
        c            = log(rels);
        goodw        = 1:nd; %1:49 to just use dates before 3rd storm
        goodp        = abs(c)>0.2;
        for i=1:nrx*2+1
            for j=1:nry*2+1
                
                phsmod       = invert_phsmat(transpose(squeeze(allinthp(:,i,j))),nd);
                phsmod(isnan(phsmod))=1;
                
                d1           = angle(phsmod);d2=d1;d3=d1;
                flips        = and(goodp',d1>pi/4);
                d2(flips)    = d2(flips)-2*pi;
                flips        = and(goodp',d1<-pi/4);
                d3(flips)    = d3(flips)+2*pi;
                      
                [f1,o1]      = fit(c(goodw),d1(goodw)','poly1','weights',1./mags(goodw).^2);
                [f2,o2]      = fit(c(goodw),d2(goodw)','poly1','weights',1./mags(goodw).^2);
                [f3,o3]      = fit(c(goodw),d3(goodw)','poly1','weights',1./mags(goodw).^2);
                
                newsynth1    = f1(c);
                newsynth2    = f2(c);
                newsynth3    = f3(c);
                newsynth1(~goodp)=0; newsynth2(~goodp)=0;newsynth3(~goodp)=0;
                
                newphsmod1   = squeeze(phsmod).*conj(exp(1j*newsynth1'));
                newphsmod2   = squeeze(phsmod).*conj(exp(1j*newsynth2'));
                newphsmod3   = squeeze(phsmod).*conj(exp(1j*newsynth3'));
                neuts        = abs(mean(newphsmod1));
                negs         = abs(mean(newphsmod2));
                poss         = abs(mean(newphsmod3));
                
                resvec       = [neuts negs poss];
                
                whichi       = find(resvec==max(resvec,[],'omitnan'));
                whichi       = whichi(end);
                output.type  = whichi;
  
                switch whichi
                    case 1
                        myslcfix(:,i,j)=newsynth1;
                        slopes(i,j)=f1.p1;
                        inters(i,j)=f1.p2;
                        adjr2(i,j)=o1.adjrsquare;
                    case 2
                        myslcfix(:,i,j)=newsynth2;
                        slopes(i,j)=f2.p1;
                        inters(i,j)=f2.p2;
                        adjr2(i,j)=o2.adjrsquare;
                    case 3
                        myslcfix(:,i,j)=newsynth3;
                        slopes(i,j)=f3.p1;
                        inters(i,j)=f3.p2;
                        adjr2(i,j)=o3.adjrsquare;
                end
            end
        end
        newslcs     = slcs.*conj(exp(1j*myslcfix));
        newints     = newslcs(intid(:,1),:,:).*conj(newslcs(intid(:,2),:,:));
        newintsf    = convn(newints,windn,'same')./wgtsn;
        newcor      = abs(newintsf);
        newintsf    = newintsf./newcor;
        newintshp   = newints./conj(newintsf);
        triplets1 = trip_nm(squeeze(newintsf(:,nrx+1,nry+1)),intid,1);
         output.newhp   = newintshp(:,nrx+1,nry+1);
        output.newlp   = newintsf(:,nrx+1,nry+1);
        %output.newphs  = newints(:,nrx+1,nry+1);
        output.newcors = newcor(:,nrx+1,nry+1);
        output.trip1=triplets1;
        output.slope=slopes(nrx+1,nry+1);
        output.inter=inters(nrx+1,nry+1);
        output.adjr2=adjr2(nrx+1,nry+1);
       if(plotflag)
            figure
            subplot(2,2,1)
            triplot(output.cors,dn);
            cax=caxis;
            subplot(2,2,2)
            triplot(output.newcors,dn);
            caxis(cax)
            subplot(2,2,3)
            triplot(angle(output.trip0),dn);
            cax=caxis;
            subplot(2,2,4)
            triplot(angle(output.trip1),dn);
            caxis(cax)
        end
    else
        output.type=0;
        output.newhp=[];
        output.newlp=[];
        %output.newphs=[];
        output.newcors=[];
        output.trip1 = [];
        output.slope=[];
        output.inter=[];
        output.adjr2=[];
    end
else
    output.raincount=0;
    output.type=0;
        output.newhp=[];
        output.newlp=[];
        %output.newphs=[];
        output.newcors=[];
        output.trip1 = [];
        output.slope=[];
        output.inter=[];
        output.adjr2=[];
end


return
if(plotflag)
    
    badi=[];
    
    figure('name',[num2str(xpt) ' ' num2str(ypt)])
    colormap('jet')
    for l=1:length(output)
        subplot(3,4,l*2-1)
        badi=[];
        dn=[output(l).dn];
        nd=length(dn);
        jnk=nan(nd);
        
        jnk([output(l).cids])=[output(l).cors];
        pcolor(dn,dn,jnk'),shading flat,set(gca,'ydir','reverse');
        hold on
        fixplot
        
        title(['cors' pols{l}])
    end
    subplot(3,1,2)
    hold on
    cols={'k','r'};
    for i=1:length(output)
        plot([output(i).dn],[output(i).mags],cols{k})
    end
    axis tight
    ax=axis;
    datetick('x','mmmYY')
    axis(ax);
    for l=1:length(output)
        dn=[output(l).dn];
        id1=[output(l).id1];
        id2=[output(l).id2];
        rels=[output(l).rels];
        modp=[output(l).modp];
        c0=[output(l).c0];
        nd=length(dn);
        ni=length(id1);
        
        Gi=zeros(ni,nd-1); %intervals, perm cor loss
        for i=1:ni
            Gi(i,id1(i):id2(i)-1)=1;
        end
        
        Gr=zeros(ni,nd); %rel cor on dates
        for i=1:ni
            Gr(i,id1(i))=1;
            Gr(i,id2(i))=-1;
        end
        
        goodrel = isfinite(rels);
        s2=exp(Gi*log(modp));
        s3=exp(-abs(Gr(:,goodrel)*log(rels(goodrel))));
        synth=c0*s2.*s3;
        subplot(3,4,l*2)
        jnk=nan(nd);
        jnk([output(l).cids])=synth;
        pcolor(dn,dn,jnk'),shading flat,set(gca,'ydir','reverse');
        hold on
        fixplot
        
        title(['synth' pols{l}])
        
        
        dni   = dn(1:nd-1)+diff(dn)/2;
        
        subplot(3,1,3)
        plot(dn,rels,'-','color',cols{l});
        hold on
        plot(dni,modp,'--','color',cols{l});
        plot(dn,c0*ones(size(dn)),':','color',cols{l});
        
    end
    axis tight
    ax=axis;
    datetick('x','mmmYY')
    axis(ax);
    
    figure
    subplot(2,3,1)
    triplot(output(1).hp,output(1).dn)
    subplot(2,3,4)
    triplot(output(2).hp,output(2).dn)
   
    subplot(2,3,2)
    for i=1:length(output)
        
        dn=[output(i).dn];
        nd=length(dn);
        jnk=nan(nd);
        jnk([output(i).cids])=exp(1j*output(i).phs);
        jnkang=diag(jnk(1:end-1,1:end-1).*jnk(2:end,2:end).*conj(jnk(2:end,1:end-1)));
        plot(dn(2:end-1),angle(jnkang(1:end-1)),'-','color',cols{i});
        hold on
        axis tight
        ax=axis;
        datetick('x','mmmYY')
        axis(ax)
        title(['shortest phase triplet closure, ' pols{i}]);
    end
    ax=axis;
    if(exist('dnr','var'))
        for i=1:length(dnr)
            plot([dnr(i) dnr(i)],[ax(3) ax(4)],'m')
        end
    end
    
    subplot(2,3,5)
    for i=1:length(output)
        
        dn=[output(i).dn];
        nd=length(dn);
        jnk=nan(nd);
        jnk([output(i).cids])=exp(1j*output(i).phs);
        jnkang=diag(jnk(1:end-1,1:end-1).*jnk(2:end,2:end).*conj(jnk(2:end,1:end-1)));
        plot(dn(2:end-1),cumsum(angle(jnkang(1:end-1))),'-','color',cols{i});
        hold on
        axis tight
        ax=axis;
        datetick('x','mmmYY')
        axis(ax)
        title(['cumulative sum shortest phase triplet closure, ' pols{i}]);
    end
    ax=axis;
    if(exist('dnr','var'))
        for i=1:length(dnr)
            plot([dnr(i) dnr(i)],[ax(3) ax(4)],'m')
        end
    end
    if(length(output)==2)
        subplot(2,3,3)
        
        dn1        = output(1).dn;
        dn2        = output(2).dn;
        [~,i1,~] = intersect(dn1,dn2);
        ints    = [output(1).id1 output(1).id2];
        gint    = sum(ismember(ints,i1),2)==2;
        phs1    = exp(1j*output(1).phs(gint));
        phs2    = exp(1j*output(2).phs);
        phsdiff = angle(phs1.*conj(phs2));
        
        triplot(phsdiff,dn2)
        hold on
        fixplot
        caxis([-pi pi])
        title('VV VH phs diff')
    end
    subplot(2,3,6)
    plot(output(1).hp,output(1).cors,'.')
end


