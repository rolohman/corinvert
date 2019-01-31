function [cors,allrels,allmodp,allc0,synth]= plot_fit(xpt,ypt,tdir)
decide_ints_stack

if(nargin==2)
    tdir='results_dates';
end

dni   = dn(1:nd-1)+diff(dn)/2;
%dni=dni(1:end-1);

intdt=diff(dn);
Gi=zeros(ni,nd-1); %intervals, perm cor loss
for i=1:ni
    Gi(i,id1(i):id2(i)-1)=1;
end
%Gi=Gi(:,1:end);

Gr=zeros(ni,nd); %rel cor on dates
for i=1:ni
    Gr(i,id1(i))=1;
    Gr(i,id2(i))=-1;
end


%matrix that maps ints to a grid
cids=[];
for i=1:nd
    cids=[cids (nd)*(i-1)+[i:nd-1]];
end

allrels=zeros(nd,length(xpt));
allmodp=zeros(nd-1,length(xpt));
allc0=zeros(1,length(xpt));

for k=1:length(xpt)
    cors=zeros(ni,1);
    rels=zeros(nd,1);
    modp=zeros(nd-1,1);
    %open data files
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
    cors=zeros(ni,1);
    for i=1:ni
        fseek(fidi(i),(newnx*(ypt(k)-1)+xpt(k)-1)*4,-1);
        cors(i)=fread(fidi(i),1,'real*4');
    end
    fclose('all');
    
    %load data and models
    fid0=fopen([tdir '/c0.cor'],'r');
    fseek(fid0,(newnx*(ypt(k)-1)+xpt(k)-1)*4,-1);
    c0=fread(fid0,1,'real*4')
    fclose(fid0);
    for i=1:nd
        fidr(i)=fopen([tdir '/rel_' dates(i).name '.cor'],'r');
    end
    rels=zeros(nd,1);
    for i=1:nd
        fseek(fidr(i),(newnx*(ypt(k)-1)+xpt(k)-1)*4,-1);
        rels(i)=fread(fidr(i),1,'real*4');
    end
    fclose('all');
    for i=1:nd-1
        fidp(i)=fopen([tdir '/perm_' dates(i).name '_' dates(i+1).name '.cor'],'r');
    end
    modp=zeros(nd-2,1);
    for i=1:nd-1
        fseek(fidp(i),(newnx*(ypt(k)-1)+xpt(k)-1)*4,-1);
        modp(i)=fread(fidp(i),1,'real*4');
    end
    modp(isnan(modp))=0;
    fclose('all');
    
    good        = cors>0;
    cors(~good) = NaN;
    ngood       = sum(good,1);
    
    badint  = isnan(cors(:));
    goodrel = isfinite(rels);
    
    s2=exp(Gi*log(modp));
    s3=exp(-abs(Gr(:,goodrel)*log(rels(goodrel))));
    synth=c0*s2.*s3;
    
    res           = cors-synth;
    synth(badint) = NaN;
    
    jnk=nan(nd);
    figure('name',[num2str(xpt(k)) ' ' num2str(ypt(k))])
    colormap('jet')
    badi=[];
    subplot(2,2,1)
    jnk(cids)=cors;
    pcolor(dn,dn,jnk'),shading flat,set(gca,'ydir','reverse');
    hold on
    fixplot
    title('data')
    
    subplot(2,2,2)
    jnk(cids)=synth;
    pcolor(dn,dn,jnk'),shading flat,set(gca,'ydir','reverse');
    hold on
    fixplot
    title(norm(cors(~badint)-synth(~badint)));
    
    subplot(2,1,2)
    l(1)=plot(dn,[rels],'b.-');
    hold on
    l(2)=plot(dni,modp,'g.-');
    l(3)=plot(dn,c0*ones(size(dn)),'c');
    axis tight
    ax=axis;
    
    %legend(l,'cr','cp','c0')
    datetick('x','mmmYY')
    axis(ax)
    title('model vs. time')
    
    allc0(k)=c0;
    allmodp(:,k)=modp;
    allrels(:,k)=rels;
    
end