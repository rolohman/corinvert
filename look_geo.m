dirs={'T54','T156','T47'};
pol='_VV';
r=4;
a=4;
dirs=dirs(3);
for i=1:length(dirs)
    d=[dirs{i} '/geo' pol '/'];
    chdir(d);
    files=dir('*.geo');
    for j=1:length(files)
       % if(~regexp(files(j).name,'4r_4a'))
        newfile=[files(j).name(1:end-8) '_' num2str(r) 'r_' num2str(a) 'a.cor.geo'];
        if(~exist(newfile,'file'))
        
            
            fidi=fopen(files(j).name,'r');
            fidr=fopen('rows.geo','r');
            fido=fopen('tmp','w');
            for k=1:ny_geo
                tmp1=fread(fidi,nx_geo,'real*4');
                tmp2=fread(fidr,nx_geo,'real*4');
                tmp1(tmp2==0)=NaN;
                fwrite(fido,tmp1,'real*4');
            end
            movefile('tmp',files(j).name);
            fclose('all');
            command=['looks.py -i ' files(j).name ' -o ' newfile ' -r ' num2str(r) ' -a ' num2str(a)];
            system(command);
        else
            disp(['already looked down ' newfile])
        end
        %end
    end
    chdir('../../');
end
