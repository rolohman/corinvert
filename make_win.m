function [windx,windy]=make_win(rx,ry,windowtype);
switch windowtype
    case 0
        windx=[1/(rx+1):1/(rx+1):1 (1-1/(rx+1)):-1/(rx+1):1/(rx+1)];
        windy=[1/(ry+1):1/(ry+1):1 (1-1/(ry+1)):-1/(ry+1):1/(ry+1)];
    case 1
        windx=zeros(1,rx*2+1); windx(floor(rx/2):ceil(rx/2)+rx)=1;
        windy=zeros(1,ry*2+1); windy(floor(ry/2):ceil(ry/2)+ry)=1;
    case 2
        windx=exp(-(-rx:rx).^2/2/(rx/2)^2);
        windy=exp(-(-ry:ry).^2/2/(ry/2)^2);
end
windx=windx/sum(windx);
windy=windy/sum(windy);