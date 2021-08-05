function grad = vertdsdtgrad(chnkr)
%VERTDSDTGRAD gradient of chunk arclength density with respect to each 
% vertex
%
% Syntax: 
%   grad = vertgrad(chnkr)
%
% Input:
%
% chnkr - chunk polygon object
%
% Output:
%
% grad - is a sparse 
%           (chnkr.npts) x (chnkr.dim \cdot chnkr.nvert) 
%        matrix which gives the gradient of the chunker coordinates with 
%        respect to the vertex locations
%


nvert = chnkr.nvert;
dim = chnkr.dim;
npt = chnkr.npt;

grad = sparse(npt,dim*nvert);

itmp1 = [0,(1:(nvert-1))];
itmp2 = 1:nvert; 

isgrad = [dim^2*nvert+1;dim^2*nvert] + dim * [itmp1;itmp2];
igall = chnkr.igraddatarows;
isg = igall(isgrad);

for i = 1:nvert
    jshift = (i-1)*chnkr.dim;
    j1 = 1:dim; j1 = j1(:); 
    
    ipts = 1:chnkr.npt; ipts = ipts(:);
    ii = repmat(ipts.',dim,1); ii = ii(:);
    n1 = numel(ipts);
    jj = repmat(j1,n1,1); jj = jshift + jj(:);
    isgi = isg(1,i):isg(2,i);
    vv = chnkr.data(isgi,:,:); vv = vv(:);
    
    grad = grad + sparse(ii,jj,vv,npt,dim*nvert);    
    
end

