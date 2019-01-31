function [difamps,allamps]= plot_amp_stack(xpt,ypt,pol)
if(pol==1)
	pol='_VV';
else
	pol='_VH';
end
pol
decide_ints_stack
%matrix that maps ints to a grid
cids=[];
for i=1:nd
    cids=[cids (nd)*(i-1)+[i:nd-1]];
end


for k=1:length(xpt)
    amps=zeros(ni,1);
    %open data files
    for i=1:nd
        slcdir=['merged/SLC' pol '/' dates(i).name '/' ];
        slcfile=[slcdir  dates(i).name '.slc.full'];
        magfile=[slcdir dates(i).name '_7_3.mag'];
        if(and(exist(slcfile,'file'),~exist(magfile,'file')))
            command=['imageMath.py -e=''abs(a)'' -o ''' slcdir dates(i).name '.mag'' -t float --a=''' slcfile ''''];
            system(command);
            command=['looks.py -i ' slcdir dates(i).name '.mag -o ' magfile ' -r 7 -a 3'];
            system(command);
           
        end
        
        if(exist(magfile,'file'))
            fidi(i)=fopen(magfile,'r');
        else
            disp(['where is ' magfile '?'])
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
