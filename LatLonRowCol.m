function [col,row,x1,y1]= LatLonRowCol(x,y,colf,rowf,latlonflag);
%row,column in downlooked SLC, plus x1,y1 pixel location in row/col file.
%Usually "rows.geo".
%latlonflag
%case 1 = x,y given in lat/lon, x1,y1 given in pixels
%case 2 = x,y given in pixels (from rowf,colf file), x1,y1, given in latlon


if(exist(rowf))
    [~,b]=system(['grep rasterXSize ' rowf '.vrt']);
    tmp=regexp(b,'rasterXSize="(\d+)" rasterYSize="(\d+)">','tokens');
    if(length(tmp)==1)
        nx=str2num(tmp{1}{1});
        ny=str2num(tmp{1}{2});
    else
        disp('nx not found')
    end
    [~,b]=system(['grep GeoTransform ' rowf '.vrt']);
    tmp=regexp(b,'[,<>]','split');
    if(length(tmp)>1)
        lon1=str2num(tmp{3});
        dlon=str2num(tmp{4});
        lat1=str2num(tmp{6});
        dlat=str2num(tmp{8});
    else
        disp('GeoTransform not found')
    end
  else
    disp([rowf ' doesn''t exist']);
end

switch latlonflag
    case 1
        disp('converting x,y in lat/lon to pixels in slc')
        %figure out which row, pixel in row,col.geo is the right lat/lon
        x1 = floor((x-lon1)/dlon);
        y1 = floor((y-lat1)/dlat);
        %read the slc row/col values
        fid=fopen(rowf,'r');
        fseek(fid,(nx*(y1-1)+x1-1)*4,-1);  %one pixel before goal pixel
        row=fread(fid,1,'real*4');
        fclose(fid);
        
    case 2        
        disp('converting pixel locations in geo file to pixels in slc');
       %figure out what the latlon vals are for a pixel in row,col.geo 
       x1 = x*dlon + lon1;
       y1 = y*dlat + lat1;  
       %think about corner locations if using for anything important.
       %read the slc row/col values
       fid=fopen(rowf,'r');
       fseek(fid,(nx*(y-1)+x-1)*4,-1);  %one pixel before goal pixel
       row=fread(fid,1,'real*4');
       fclose(fid);
       fid=fopen(colf,'r');
       fseek(fid,(nx*(y-1)+x-1)*4,-1);  %one pixel before goal pixel
       col=fread(fid,1,'real*4');
       fclose(fid);
end
        