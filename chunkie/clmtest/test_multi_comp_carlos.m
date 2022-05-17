% This file is supposed to be a multiple component test for solving
% the transmission problem in a multiply connected geometry


centers = [-0.5,-0.5,0.3; 0.5,-0.5,0];
phis = [0;0.0;0];
scales = [0.3;0.3;0.3];
[~,ncurve] = size(centers);

kin = 3.0;
kout = 5.0;



% set target locations
%receptors
r_tgt = 10;
n_tgt = 100;
t_tgt = 0:2*pi/n_tgt:2*pi-2*pi/n_tgt;
x_t   = r_tgt * cos(t_tgt);
y_t   = r_tgt * sin(t_tgt);    
tgt   = [ x_t; y_t];

% Set interior point
src0 = [0.01;-1.1];%horseshoe


% incidence directions
n_dir = 88;
t_dir = 0:2*pi/n_dir:2*pi-2*pi/n_dir;
x_dir = cos(t_dir);
y_dir = sin(t_dir);
dir =[ x_dir; y_dir ];



ppw = 20;
% interior point in unscaled curve
src0 = [-0.1;0.12];

pref = []; 
pref.k = 16;


cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 0;

narms = 3;
amp = 0.25;
start = tic; chnkr0 = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref);
wts = weights(chnkr0);
L0 = sum(wts(:));
t1 = toc(start);

chnkr(1,ncurve) = chunker();


ndomain = ncurve+1;
src = zeros(2,ndomain);



k = zeros(ndomain,1);
k(1:ncurve) = kin;
k(ncurve+1) = kout;

kmax = max(k);

nch0 = max(10,ceil(ppw*L0*kmax/2/pi/pref.k));



% TE mode
coef = ones(size(k));

%TM mode
%coef = kout^2/kin^2;

c = zeros(2,ncurve);
c(1,:) = 1:ncurve;
c(2,:) = ndomain;



k1 = zeros(1,ncurve); % wave numbers for the interior domain
k2 = zeros(1,ncurve); % wave numbers for the exterior domain

for i=1:ncurve
  k1(i) = k(c(1,i));
  k2(i) = k(c(2,i));
end

fcurve = cell(1,ncurve);
cparams = cell(1,ncurve);

for icurve=1:ncurve
    fcurve{icurve} = @(t) starfish(t,narms,amp,centers(1:2,icurve), ...
      phis(icurve),scales(icurve));
    src(:,icurve) = src0*scales(icurve) + centers(:,icurve);
    cparams{icurve}.ta = 0;
    cparams{icurve}.tb = 2*pi;
    cparams{icurve}.ifclosed = true;
end

src(1,ndomain) = 10.1;
src(2,ndomain) = 3.1;

nch = ones(ncurve,1)*nch0;
alpha1 = zeros(1,ncurve);
alpha2 = zeros(1,ncurve);
for icurve=1:ncurve
  chnkr(icurve) = chunkerfuncuni(fcurve{icurve},nch(icurve),cparams{icurve},pref);
  chnkr(icurve) = chnkr(icurve).makedatarows(4);
  chnkr(icurve).data(1,:,:) = k(c(1,icurve));
  chnkr(icurve).data(2,:,:) = k(c(2,icurve));
  chnkr(icurve).data(3,:,:) = coef(c(1,icurve));
  chnkr(icurve).data(4,:,:) = coef(c(1,icurve));
  alpha1(icurve) = 2/(coef(c(1,icurve))+coef(c(2,icurve)));
  alpha2(icurve) = 2/(1/coef(c(1,icurve))+1/coef(c(2,icurve)));
end
ngl = chnkr(1).k;
np = sum(nch(1:ncurve))*ngl;

% form system matrix

% number of Gauss-Legendre nodes on each chunk
ngl = 16;
[glnodes,glwts] = lege.exps(ngl);


% Treat the representation as if it were a 2x2 operator so that four layer 
% potentials D, S, D', S' can be evaluated together. This will avoid 
% redundant evaluation of Hankel functions.
opdims(1)=2;opdims(2)=2;

% set up GGQ machinery
[logquad] = chnk.quadggq.setuplogquad(ngl,opdims);


% build the system matrix
rpars = [];
rpars.k = k;
rpars.c = c;
rpars.coef = coef;

disp(' ')
disp('Step 1: build the system matrix directly.')
start = tic;
ilist=[];
opts = [];
[M,np,~,~] = clm.buildmat_fast(chnkr,rpars,opts,opdims,glwts,ilist,logquad);
M = M + eye(2*np);
dt = toc(start);

disp(['System matrix construction time = ', num2str(dt), ' seconds'])





% construct artificial boundary data for testing purpose

rhs = zeros(2*np,1);

rhsinc = zeros(2*np,n_dir);
chnkrtotal = merge(chnkr);

for i=1:ncurve
  d1 = c(1,i); % interior domain index
  d2 = c(2,i); % exterior domain index

  j1 = d1 + 1; % src index for the interior domain
  if j1 > ndomain
    j1 = j1 - ndomain;
  end
  j2 = d2 + 1; % src index for the exterior domain
  if j2 > ndomain
    j2 = j2 - ndomain;
  end
  
  fprintf('icurve: %d     d1: %d     d2 :%d\n',i,d1,d2);

  c1 = coef(d1);
  c2 = coef(d2);

  ind1 = sum(nch(1:i-1))*ngl*2+(1:2:2*nch(i)*ngl);
  ind2 = sum(nch(1:i-1))*ngl*2+(2:2:2*nch(i)*ngl);

  targnorm = chnkr(i).n;
  nx = targnorm(1,:); nx=nx.';
  ny = targnorm(2,:); ny=ny.';

  [u1,gradu1]=chnk.helm2d.green(k1(i),src(:,j1),chnkr(i).r(:,:));
  du1dn = gradu1(:,1).*nx + gradu1(:,2).*ny;

  [u2,gradu2]=chnk.helm2d.green(k2(i),src(:,j2),chnkr(i).r(:,:));
  du2dn = gradu2(:,1).*nx + gradu2(:,2).*ny;

  rhs(ind1) = -alpha1(i)*(u1-u2);
  rhs(ind2) =  alpha2(i)*(1/c1*du1dn-1/c2*du2dn);
  xs = chnkr(i).r(1,:);
  ys = chnkr(i).r(2,:);
  dxs = chnkr(i).d(1,:);
  dys = chnkr(i).d(2,:);
  ds = sqrt(dxs.^2 + dys.^2);
  
  nn = length(xs);
  
  
  
  u1inc = zeros(nn,n_dir);
  u2inc = zeros(nn,n_dir);
  dudn1inc = zeros(nn,n_dir);
  dudn2inc = zeros(nn,n_dir);
  if(d1 == ndomain)
    u1inc  = exp(1i *k1(i) * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    dudn1inc = 1i* k1(i) * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            u1inc;
  end
  if(d2 == ndomain)
      u2inc  = exp(1i *k2(i) * (bsxfun(@times,xs',x_dir)+bsxfun(@times,ys',y_dir)));
    dudn2inc = 1i* k2(i) * (bsxfun(@times,dys',x_dir)-bsxfun(@times,dxs',y_dir))./repmat(ds',1,n_dir) .* ...
            u2inc;
      
  end
  rhsinc(ind1,:) = -alpha1(i)*(u1inc-u2inc);
  rhsinc(ind2,:) =  alpha2(i)*(1/c1*dudn1inc-1/c2*dudn2inc);

end


Minv = inv(M);
sol = Minv*rhs;
targ = src;
%sol1 = sol(1:2:2*np); % double layer density
%sol2 = sol(2:2:2*np); % single layer density
%abs(sol(1))


uexact = zeros(ndomain,1);
ucomp = zeros(ndomain,1);

% compute the exact solution
for i=1:ndomain
  j=i+1;
  if j > ndomain
    j = j - ndomain;
  end
  
  uexact(i) = uexact(i) + chnk.helm2d.green(k(i),src(:,j),targ(:,i));
end

% compute the numerical solution
for i=1:ndomain
  evalkern = @(s,t) chnk.helm2d.kern(k(i),s,t,'eval',coef(i));
  ucomp(i) = chunkerkerneval(chnkrtotal,evalkern,sol,targ(:,i));
end

uerror = abs(ucomp-uexact)./abs(uexact);
disp(' ')
disp('Now check the accuracy of numerical solutions')
disp('Exact value               Numerical value           Error')  
fprintf('%0.15e     %0.15e     %7.1e\n', [real(uexact).'; real(ucomp).'; real(uerror)'])

sol_pw = Minv*rhsinc;
srcinfo = [];
srcinfo.r = chnkrtotal.r;
srcinfo.d = chnkrtotal.d;
srcinfo.n = chnkrtotal.n;
srcinfo.r = reshape(srcinfo.r,2,chnkrtotal.k*chnkrtotal.nch);
srcinfo.d = reshape(srcinfo.d,2,chnkrtotal.k*chnkrtotal.nch);
srcinfo.n = reshape(srcinfo.n,2,chnkrtotal.k*chnkrtotal.nch);

targinfo = [];
targinfo.r = tgt;

xkern = chnk.helm2d.kern(k(ndomain),srcinfo,targinfo,'eval',coef(ndomain));
wts = weights(chnkrtotal);
wts = reshape(wts,chnkrtotal.k*chnkrtotal.nch,1);
wts_rep = repmat(wts,[1,n_dir]);
sol_pw(1:2:end,:) = sol_pw(1:2:end,:).*wts_rep;
sol_pw(2:2:end,:) = sol_pw(2:2:end,:).*wts_rep;
umeas = xkern*sol_pw;


