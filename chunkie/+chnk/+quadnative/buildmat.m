function submat = buildmat(chnkr,kern,opdims,i,j,wts)
%CHNK.QUADSMOOTH.BUILDMAT build matrix for far interactions with this kernel
% assuming that the smooth rule is sufficient
% 
%


if nargin < 4
    i = 1:chnkr.nch;
end
if nargin < 5
    j = 1:chnkr.nch;
end
if nargin < 6
    [~,wts] = lege.exps(chnkr.k);
end

% grab specific boundary data

r = chnkr.r;
d = chnkr.d;
n = chnkr.n;
d2 = chnkr.d2;
h = chnkr.h;
if(chnkr.hasdata)
    data = chnkr.data;
end

[dim,k,~] = size(r);
rs = r(:,:,j); rt = r(:,:,i); 
ds = d(:,:,j); dt = d(:,:,i); 
ns = n(:,:,j); nt = n(:,:,i);
d2s = d2(:,:,j); d2t = d2(:,:,i); 
if(chnkr.hasdata)
    dds = data(:,:,j); ddt = data(:,:,i);
else
    dds = []; ddt = [];
end

rs = reshape(rs,dim,k*length(j)); rt = reshape(rt,dim,k*length(i));
ds = reshape(ds,dim,k*length(j)); dt = reshape(dt,dim,k*length(i));
ns = reshape(ns,dim,k*length(j)); nt = reshape(nt,dim,k*length(i));
d2s = reshape(d2s,dim,k*length(j)); d2t = reshape(d2t,dim,k*length(i));
if(chnkr.hasdata)
    [ddim,~,~] = size(data);
    dds = reshape(dds,ddim,k*length(j)); ddt = reshape(ddt,ddim,k*length(i));
end


srcinfo = []; srcinfo.r = rs; srcinfo.d = ds; srcinfo.d2 = d2s; srcinfo.n = ns; srcinfo.data = dds;
targinfo = []; targinfo.r = rt; targinfo.d = dt; targinfo.d2 = d2t; targinfo.n = nt; targinfo.data = ddt;
hs = h(j); ht = h(i);

dsnrms = sqrt(sum(ds.^2,1));
%dsnrms = ds(1,:,:); % for complex contour, by SJ 09/30/21
%taus = bsxfun(@rdivide,ds,dsnrms);

%dtnrms = sqrt(sum(dt.^2,1));
%taut = bsxfun(@rdivide,dt,dtnrms);

ws = kron(hs(:),wts(:));

dsdt = dsnrms(:).*ws;

dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
dsdtndim2 = dsdtndim2(:);

submat = bsxfun(@times,kern(srcinfo,targinfo),(dsdtndim2(:)).');

end
