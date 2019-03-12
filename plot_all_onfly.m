function [output]= plot_all_onfly(xpt,ypt,plotflag,slcflag)
if(~exist('slcflag','var'))
    slcflag=0;
end
if(~exist('plotflag','var'))
    plotflag=0;
end
      
lonfile='merged/geom_master/lon.rdr.4alks_15rlks.full';
latfile='merged/geom_master/lat.rdr.4alks_15rlks.full';



pols={'_VV','_VH'};
for l=1:length(pols)
    pol=pols{l}
    clear dates
    decide_ints_stack
    if(nd==0)
        disp(['no dates for ' pol])
    else
        output(l,1).dn=dn;
        output(l,1).id1=id1;
        output(l,1).id2=id2;
                 
        rx=rlooks;
        ry=alooks;
        
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
        slc1 = zeros(ry*2+1,nx);
        slc2 = slc1;
        
        
        dni   = dn(1:nd-1)+diff(dn)/2;
        %dni=dni(1:end-1);
        
        intdt=diff(dn);
        
        
        %matrix that maps ints to a grid
        cids=[];
        for i=1:nd
            cids=[cids (nd)*(i-1)+[i:nd-1]];
        end
        output(l).cids=cids;
        %open all slcs
        for i=1:nd
            dates(i).name    = files(i).name(1:8);
            dates(i).slc     = [slcdir dates(i).name '/' dates(i).name '.slc.full'];
            fidi(i)          = fopen(dates(i).slc,'r');
        end
        %check if inversion is done (assume c0 is typical)
        tdir=['results_dates' pol '/'];
        c0filename=[tdir 'c0.cor'];
        d=dir(c0filename);
        if(isempty(d))
            disp(['no inversion found for ' pol]);
            output(l).c0=NaN;
            output(l).rels=zeros(nd,1);
            output(l).modp=zeros(nd-1,1);
        else
            lines=d.bytes/newnx/4;
            disp(['processed inversion up to line ' num2str(lines)]);
            for k=1:length(xpt)
                if(exist(lonfile,'file'))
                    fid=fopen(lonfile,'r');
                    fseek(fid,(newnx*(ypt(k)-1)+xpt(k)-1)*8,-1);
                    output(l,k).lon=fread(fid,1,'real*8');
                    fid=fopen(latfile,'r');
                    fseek(fid,(newnx*(ypt(k)-1)+xpt(k)-1)*8,-1);
                    output(l,k).lat=fread(fid,1,'real*8');
                end
                output(l,k).x=xpt(k);
                output(l,k).y=ypt(k);

                rels=zeros(nd,1);
                modp=zeros(nd-1,1);
                if(lines>=ypt(k))
                    fid0=fopen([tdir '/c0.cor'],'r');
                    fseek(fid0,(newnx*(ypt(k)-1)+xpt(k)-1)*4,-1);
                    output(l,k).c0=fread(fid0,1,'real*4');
                    fclose(fid0);
                    for i=1:nd
                        fidr(i)=fopen([tdir '/rel_' dates(i).name '.cor'],'r');
                    end
                    for i=1:nd
                        if(fidr(i)>0)
                            fseek(fidr(i),(newnx*(ypt(k)-1)+xpt(k)-1)*4,-1);
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
                            fseek(fidp(i),(newnx*(ypt(k)-1)+xpt(k)-1)*4,-1);
                            modp(i)=fread(fidp(i),1,'real*4');
                            fclose(fidp(i));
                        else
                            modp(i)=1;
                        end
                        
                    end
                    modp(isnan(modp))=1;
                    
                else
                    disp(['no model yet for line ' num2str(ypt(k))]);
                    output(l,k).c0=NaN;
                end
                
                output(l,k).rels=rels;
                output(l,k).modp=modp;
            end
        end
        for k=1:length(xpt)
            clear cors phs
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
            mags=reshape(slcs,nd,(rx*2+1)*(ry*2+1));
            mags=abs(mags);
            mags=median(mags,2,'omitnan');
            output(l,k).mags=mags;
            if(slcflag)
                output(l,k).slcs=slcs;
            end
            for i=1:ni
                slc1=squeeze(slcs(id1(i),:,:));
                slc2=squeeze(slcs(id2(i),:,:));
              
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
                cpx3 = cpx3(rx+1);
  
                sm   = abs(cpx3);
                sm(isnan(sm))=0;
                
                phs(i)=angle(cpx3);
                cors(i)=sm;
            end
            
            good        = cors>0;
            cors(~good) = NaN;
            ngood       = sum(good,1);
            
            
            output(l,k).phs=phs;
            output(l,k).cors=cors;
            
        end
        fclose('all');
    end
end

if(plotflag)
    
    badi=[];
    
    figure('name',[num2str(xpt(k)) ' ' num2str(ypt(k))])
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
    subplot(1,2,1)
    for i=1:length(output)
        
        dn=[output(i).dn];
        nd=length(dn);
        jnk=nan(nd);
        jnk([output(i).cids])=exp(im*output(i).phs);
        jnkang=diag(jnk(1:end-1,1:end-1).*jnk(2:end,2:end).*conj(jnk(2:end,1:end-1)));
        plot(dn(2:end-1),angle(jnkang(1:end-1)),'-','color',cols{i});
        hold on
        datetick('x','mmmYY')
        title(['shortest phase triplet closure, ' pols{i}]);
    end
    if(length(output)==2)
        subplot(1,2,2)
        
        dn1        = output(1).dn;
        dn2        = output(2).dn;
        [~,i1,~] = intersect(dn1,dn2);
        ints    = [output(1).id1 output(1).id2];
        gint    = sum(ismember(ints,i1),2)==2;
        phs1    = exp(im*output(1).phs(gint));
        phs2    = exp(im*output(2).phs);
        phsdiff = angle(phs1.*conj(phs2));
        
        triplot(phsdiff,dn2)
        fixplot
        caxis([-pi pi])
        title('VV VH phs diff')
    end
end


