decide_ints_stack

adir     = 'results_dates/';
if(~exist(adir,'dir'))
    disp(['creating directory ' adir]);
    mkdir(adir)
end
count=0;
for i=1:ni
    cordir=(['cordir/' dates(id1(i)).name '/']);
    intcordir=[cordir dates(id2(i)).name '/'];
    corfile_small=[cordir dates(id1(i)).name '_' dates(id2(i)).name '_' num2str(rlooks) 'rlk_' num2str(alooks) 'alk.cor'];
    if(exist(corfile_small,'file'))
       count=count+1;
        fidi(count)=fopen(corfile_small,'r');
    else
%        disp(['where is ' corfile_small '?'])
    end
end
disp([num2str(count) ' premade cor files']);
ni=count;
alpha=50;

fid1=fopen([adir 'cmax_robust.cor'],'w');
fid2=fopen([adir 'cmax.cor'],'w');
fid3=fopen([adir 'cmin_robust.cor'],'w');
fid4=fopen([adir 'cmin.cor'],'w');

for j=1:newny
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


