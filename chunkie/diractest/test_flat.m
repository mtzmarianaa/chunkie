% This function tests the flat interface problem
% Representation S [\tilde{\rho}]
%

m = 2.0;
E = 0.01;

close('all')


xmin = -60;
xmax = 60;
cparams = [];
cparams.ta = xmin;
cparams.tb = xmax;
cparams.ifclosed = 0;

nch = 100;
chnkr = chunkerfuncuni(@(t) flatinterface(t),...
   nch,cparams);


zk = 1j*sqrt(m^2-E^2);
fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'s',1);
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmat1 = chunkermat(chnkr,fkern,opts);
t1 = toc(start);


zk = E;
fkern2 = @(s,t) chnk.helm1d.kern(zk,s,t,'s');
opdims(1) = 1; opdims(2) = 1;
opts = [];
start = tic; sysmat2 = chunkermat(chnkr,fkern2,opts)/E;
t1 = toc(start);


n = chnkr.k*chnkr.nch;
mat1 = eye(n)-2*m*sysmat1;
mat2 = eye(n)+sysmat2;

[~,s1,~] = svd(mat1);
[~,s2,~] = svd(mat2);


xs = chnkr.r(1,:);
xs = xs(:);
xx = xs(:);
ys = chnkr.r(2,:);
ys = ys(:);
yy = ys(:);
kh = 1j*sqrt(m^2-E^2);

x0 = 0;
y0 = 1.0;
source = [x0 y0];
rr = sqrt((xs-source(1)).^2+(ys-source(2)).^2);
rr = rr(:);
u_test = (1i/4)*besselh(0,1,kh*rr);






mat = mat1*mat2;

ich1 = 10;
ich2 = nch-10;
istart = (ich1-1)*chnkr.k +1;
iend = (ich2-1)*chnkr.k;

matt = mat(istart:iend,istart:iend);
uu = -2*m*u_test(istart:iend);

sol = matt\uu;

[~,s3,~] = svd(matt);




figure()
semilogy(diag(s1),'k.')
ylim([10^-6,10])
figure()
semilogy(diag(s2),'k.')
ylim([10^-6,10])
figure()
semilogy(diag(s3),'k.')
ylim([10^-6,10])

tt = chnkr.data;

figure()
plot(tt(istart:iend),real(sol),'k.')


sol2 = zeros(n,1);
sol2(istart:iend) = sol;
sol3 = mat2*sol2;
figure()
plot(tt(:),real(sol3),'k.')

sol4 = sysmat1*sol3 - u_test;

figure()
plot(chnkr)

return

vuse = conj(vs(:));
vvals = vuse./sol4;
vvals2 = vuse-sol4;
figure(10)
clf()
plot(tt(istart:iend),abs(vvals(istart:iend)),'k.')

figure(11)
clf()
plot(tt(istart:iend),abs(vvals2(istart:iend)),'k.')
