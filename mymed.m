function value=mymed(crs,wgts)

[nd,nx]    = size(crs);

if(nd==1)
    value=crs;
else
    [s,sortid] = sort(crs);
    
    if(nargin==2)
        ws         = wgts./repmat(sum(wgts,1,'omitnan'),nd,1);
        is         = repmat(1:nx,nd,1);
        ind        = sub2ind([nd nx],sortid(:),is(:));
        wsa        = reshape(ws(ind),[nd nx]);
        wss        = cumsum(wsa);
        [~,id]     = min(abs(wss-0.5),[],1,'omitnan');
        
        ind        = sub2ind([nd nx],id,1:nx);
        value      = s(ind);
    else
        nb=sum(isfinite(crs));
        nb=ceil(nb/2);
        gi=find(nb>0);
        ind=sub2ind([nd nx],nb(gi),gi);
        values=zeros(1,nx);
        value(gi)=s(ind);
    end
end