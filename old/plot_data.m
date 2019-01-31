function [allcors]= plot_data(xpt,ypt,plotflag)
decide_ints


dni   = dn(1:nd-1)+diff(dn)/2;
%dni=dni(1:end-1);

intdt=diff(dn);


%matrix that maps ints to a grid
cids=[];
for i=1:nd
    cids=[cids (nd)*(i-1)+[i:nd-1]];
end


for k=1:length(xpt)
    cors=zeros(ni,1);
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
    
    good        = cors>0;
    cors(~good) = NaN;
    ngood       = sum(good,1);
    allcors(:,k)=cors;
    
    if(plotflag)
badi=[];    
    jnk=nan(nd);
    figure('name',[num2str(xpt(k)) ' ' num2str(ypt(k))])
    colormap('jet')
    jnk(cids)=cors;
    pcolor(dn,dn,jnk'),shading flat,set(gca,'ydir','reverse');
hold on
fixplot

    title('data')
   end 
end