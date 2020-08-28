function [fidi,fido,online]=open_files(fidi,fido,nx,ny)
%rflag = 1, open for writing/appending, 2, only for reading
bytes=4; %need to fix for complex file inputs

%%%test output files, count lines, determine correct permissions
names    = fieldnames(fido);
testfile = getfield(fido,{1},names{1});
testfile = testfile.name;
if(exist(testfile,'file'))
    info   = dir(testfile);
    npix   = info.bytes/4; %all output files are r4 for now.  Need to /8 if slcs.
    online = floor(npix/nx);
    if(online<1)
        disp([ testfile  ' empty? starting to rewrite over all files'])
        online=0;
    elseif(online==ny)
        disp([ testfile ' already full size - move old directory to new name or delete'])
        return
    end
else
    disp('output files do not exist, starting from scratch');
    online=0;
end

if(online==0)
    filestyle='w';
else
    filestyle='a';
end

%%%Open Infiles
names=fieldnames(fidi);
nf1=length(names);
for i=1:nf1
    tmp=getfield(fidi,{1},names{i});
    nf2=length(tmp);
    for j=1:nf2
        fidt=fopen(tmp(j).name,'r'); %files to read
        fseek(fidt,online*nx*bytes,-1);
        fidi=setfield(fidi,{1},names{i},{j},'fid',fidt);
    end
end


%%% Open Outfiles
names=fieldnames(fido);
nf1=length(names);
for i=1:nf1
    tmp=getfield(fido,{1},names{i});
    nf2=length(tmp);
    for j=1:nf2
        fidt=fopen(tmp(j).name,filestyle); %files to read
        fido=setfield(fido,{1},names{i},{j},'fid',fidt);
    end
end

disp(['on line ' num2str(online)])
