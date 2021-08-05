%CHUNKERPOLYTEST

clearvars; close all;
addpaths_loc();

% pre-defined vertices for a barbell shape

verts = chnk.demo.barbell(2.0,2.0,1.0,1.0);

% for barbell(2,2,1,1) area and length are

barb_area = 9;
barb_length = 16;

% rounded corner version

cparams = [];
cparams.widths = 0.1*ones(size(verts,2),1);
cparams.eps = 1e-8;

nv = size(verts,2);
edgevals = rand(3,nv); % constant values on edges to smooth out in arclength

p.k = 16; p.dim = 2;
chnkr = chunkerpoly(verts,cparams,p,edgevals);
chnkr = chnkr.sort();
assert(checkadjinfo(chnkr) == 0);

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal

figure(2)
chnkr_ref = refine(chnkr);
clf
plot(chnkr_ref,'-x')
hold on
quiver(chnkr_ref)
axis equal

figure(3)
clf
nchplot = 1:chnkr.nch;
x = chnkr.r(1,:,nchplot);
y =chnkr.r(2,:,nchplot);
z = chnkr.data(1,:,nchplot);
plot3(x(:),y(:),z(:))

% adaptive refinement in corners (no smoothing of edge data)

p.k = 16; p.dim = 2;
cparams = [];
cparams.rounded = false;
cparams.depth = 8;
chnkr2 = chunkerpoly(verts,cparams,p);
chnkr2 = chnkr2.sort();
assert(checkadjinfo(chnkr2) == 0);

x = chnkr2.r(1,:); y = chnkr2.r(2,:); z = chnkr2.data(1,:);

figure(4)
clf
plot(chnkr2,'b-x')
hold on
quiver(chnkr2)

figure(5)
clf
plot3(x(:),y(:),z(:))

barb_area_2 = area(chnkr2);
err_area = abs(barb_area-barb_area_2)/abs(barb_area);
barb_length_2 = sum(sum(weights(chnkr2)));
err_length = abs(barb_length -barb_length_2)/abs(barb_length);
fprintf('%5.2e : diff between true/computed area\n',err_area);
fprintf('%5.2e : diff between true/computed length\n',err_length);

%

verts = randn(2,5);

cparams = [];
cparams.widths = 0.1*ones(size(verts,2),1);
cparams.autowidths = true;
cparams.autowidthsfac = 0.1;
cparams.ifclosed = 0;
cparams.eps = 1e-3;

p.k = 16; p.dim = 2;
chnkr3 = chunkerpoly(verts,cparams,p);

[ichlefts,ichrights] = incident_edges(chnkr3);

figure(6)

k = chnkr3.k;

clf
plot(chnkr3,'-r')
axis equal
hold on
nv = size(verts,2);
ivert = randi(nv);
plot(verts(1,ivert),verts(2,ivert),'go')
ichl = ichlefts{ivert};
ichr = ichrights{ivert};
ilpts =  ( (1:k).') + ( (ichl(:)-1)*k ).'; ilpts = ilpts(:);
irpts =  ( (1:k).') + ( (ichr(:)-1)*k ).'; irpts = irpts(:);
plot(chnkr3.r(1,ilpts),chnkr3.r(2,ilpts),'bx')
plot(chnkr3.r(1,irpts),chnkr3.r(2,irpts),'mx')
axis equal

%

figure(7)
chnkr_ref3 = refine(chnkr3);
clf
plot(chnkr_ref3,'-x')
hold on
quiver(chnkr_ref3)
axis equal

%


% rounded corner version


verts = chnk.demo.barbell(2.0,2.0,1.0,1.0);
cparams = [];
cparams.widths = 0.1*ones(size(verts,2),1);
cparams.eps = 1e-8;

nv = size(verts,2);

p.k = 16; p.dim = 2;
chnkr = chunkerpoly(verts,cparams,p);
chnkr = chnkr.sort();
chnkr=refine(chnkr);

% sources

ns = 10;

sources = rand(2,ns) + [3.0;3.0];

strengths = randn(ns,1);

% targets

nt = 3;

targets = 0.1*(-0.5+rand(2,nt));

kerns = @(s,t) chnk.lap2d.kern(s,t,'s');

% eval u on bdry

srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = chnkr.r;
kernmats = kerns(srcinfo,targinfo);
ubdry = kernmats*strengths;

% eval u at targets

targinfo = []; targinfo.r = targets;
kernmatstarg = kerns(srcinfo,targinfo);
utarg = kernmatstarg*strengths;

figure(8)
clf
plot(chnkr,'r')
hold on
plot(sources(1,:),sources(2,:),'bo')
plot(targets(1,:),targets(2,:),'gx')
%

% build laplace dirichlet matrix

fkern = @(s,t) chnk.lap2d.kern(s,t,'D');
start = tic; D = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(chnkr.k*chnkr.nch) + D;

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

start = tic; sol2 = sys\rhs; t1 = toc(start);

fprintf('%5.2e s : time for dense backslash solve\n',t1)

err = norm(sol-sol2,'fro')/norm(sol2,'fro');

fprintf('difference between direct and iterative %5.2e\n',err)

% evaluate at targets and compare

opts.usesmooth=false;
opts.verb=false;
opts.quadkgparams = {'RelTol',1e-16,'AbsTol',1.0e-16};
start=tic; Dsol = chunkerkerneval(chnkr,fkern,sol2,targets,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n',t1)

%

wchnkr = weights(chnkr);

relerr = norm(utarg-Dsol,'fro')/(sqrt(chnkr.nch)*norm(utarg,'fro'));
relerr2 = norm(utarg-Dsol,'inf')/dot(abs(sol(:)),wchnkr(:));
fprintf('relative frobenius error %5.2e\n',relerr);
fprintf('relative l_inf/l_1 error %5.2e\n',relerr2);

assert(relerr < 1e-10,'low precision in chunkmat test on polygon at target');
