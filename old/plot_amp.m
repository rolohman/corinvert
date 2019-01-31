function [difamps,allamps]= plot_amp(xpt,ypt)
decide_ints

%matrix that maps ints to a grid
cids=[];
for i=1:nd
    cids=[cids (nd)*(i-1)+[i:nd-1]];
end


for k=1:length(xpt)
    amps=zeros(ni,1);
    %open data files
    for i=1:nd
        slcdir=(['slcs/']);
        if(i==1)
            slcfile=[slcdir 'master_7_3.mag.full'];
        else
        slcfile=[slcdir dates(i).name '_7_3.mag.full'];
        end
        if(exist(slcfile,'file'))
            fidi(i)=fopen(slcfile,'r');
        else
            disp(['where is ' slcfile '?'])
        end
    end
    amps=zeros(nd,1);
    for i=1:nd
        fseek(fidi(i),(newnx*(ypt(k)-1)+xpt(k)-1)*4,-1);
        amps(i)=fread(fidi(i),1,'real*4');
    end
    fclose('all');
    
    difamp=amps(id1)-amps(id2);
%     jnk=nan(nd);
%     figure('name',[num2str(xpt(k)) ' ' num2str(ypt(k))])
%     colormap('jet')
%     badi=[];
%     
%     jnk(cids)=difamp;
%     pcolor(dn,dn,jnk'),shading flat,set(gca,'ydir','reverse');
%     hold on
%     %fixplot
%     title('data')
   allamps(k,:)=amps; 
difamps(k,:)=difamp;    
end