decide_ints

rx=7;
ry=3;

files=dir('2*xml');
ndates=length(files);

for i=1:length(files)
    dates(i).name    = files(i).name(1:8);    
    intdir           = [dates(testid).name '/int_' dates(testid).name '_' dates(i).name '/merged/'];
    dates(i).slc     = [intdir dates(i).name '.slc.full'];
    dates(i).origslc = [intdir 'slave.slc.full'];
end
intdir                = [dates(testid).name '/int_' dates(testid).name '_' dates(2).name '/merged/'];
dates(testid).slc     = [intdir 'master.slc.full'];
dates(testid).origslc = [intdir 'master.slc.full'];
if(~exist('cordir2'))
    mkdir('cordir2')
end


for i=1:ndates-1
    %for i=18
    for j=i+1:ndates
        cordir=(['cordir2/' dates(i).name '/']);
        if(~exist(cordir))
            mkdir(cordir)
        end
        intcordir=[cordir dates(j).name '/'];
        if(~exist(intcordir))
            mkdir(intcordir);
        end
        corfile_small=[intcordir dates(i).name '_' dates(j).name '_' num2str(rx) 'rlk_' num2str(ry) 'alk.cor'];
        if(~exist(corfile_small))
            disp(['running ' corfile_small])
            
            mycor_slcs_downlook(dates(i).slc,dates(j).slc,corfile_small,nx,ny,rx,ry,0)
        else
            disp([corfile_small ' already made'])
        end
    end
end



%end
% for i=1
%     infile=files(i).name;
%     corfile=regexprep(infile,'.int','.cor');
%     outfile=[infile '.cor'];
%     %movefile(outfile,[outfile '.old']);
%     mycor4(infile,[],outfile,nx,ny,rx,ry,0,0);
%     copyfile([corfile '.xml'],[outfile '.xml']);
%     copyfile([corfile '.vrt'],[outfile '.vrt']);
%     find_and_replace([outfile '.xml'],'cor','int.cor')
%     find_and_replace([outfile '.vrt'],'cor','int.cor')
% end


