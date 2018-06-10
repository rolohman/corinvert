


baseline_all=[];
baseline_dates=[0];
for i=1:length(ints)
  
   % system(['cd int_',num2str(ints(i,1)),'_',num2str(ints(i,2))])
file=['int_',num2str(ints(i,1)),'_',num2str(ints(i,2)),'/isce.log'];
system(['rm topophase.flat.full.xml.txt'])
system(['grep Bperp ',file,' >> topophase.flat.full.xml.txt'])

%system('sed -i -s "s/isce.log:baseline/g" topophase.flat.full.xml.txt')

 teext=textread('topophase.flat.full.xml.txt','%s');
 baseline=(str2num(char(teext(9)))+str2num(char(teext(18))))/2;
 
 baseline_all=[baseline_all; baseline];
 ints(i,3)=10^6*baseline;
 %baseline_dates=[baseline_dates;baseline_dates(i)+ baseline];
   
end


