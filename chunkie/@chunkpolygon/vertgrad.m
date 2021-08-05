function [grad] = vertgrad(chnkr)
%VERTGRAD gradient of chunk polygon nodes with respect to each vertex
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
%           (chnkr.dim \cdot chnkr.npts) x (chnkr.dim \cdot chnkr.nvert) 
%        matrix which gives the gradient of the chunker coordinates with 
%        respect to the vertex locations
%

nvert = chnkr.nvert;
dim = chnkr.dim;
npt = chnkr.npt;

grad = sparse(dim*npt,dim*nvert);

itmp1 = [0,(1:(nvert-1))];
itmp2 = 1:nvert; 

igrad = [1;0] +dim^2* [itmp1;itmp2];
isgrad = [dim^2*nvert+1;dim^2*nvert] + dim * [itmp1;itmp2];
igall = chnkr.igraddatarows;
ig = igall(igrad); ig = reshape(ig,[],nvert);
isg = igall(isgrad);

for i = 1:nvert
    jshift = (i-1)*chnkr.dim;
    j1 = 1:dim; j1 = j1(:); 
    
    ipts = 1:chnkr.npt; ipts = ipts(:);
    i1 = repmat(j1.',dim,1); i1 = i1(:);
    ii = i1 + (ipts.' - 1)*dim; ii = ii(:);
    n1 = numel(ipts);
    jj = repmat(j1,n1*dim,1); jj = jshift + jj(:);
    igi = ig(1,i):ig(2,i); igi=igi(:);
    vv = chnkr.data(igi,:,:); vv = vv(:);
    
    grad = grad + sparse(ii,jj,vv,dim*npt,dim*nvert);
    
end

