dm = 2;
de = 1;


domeg = sqrt(dm*dm-de*de);


xi = 0.51;
xis = -20:0.1:20;
vals = zeros(size(xis));

if(true)
for i=1:numel(xis)
xi = xis(i);    
zk = sqrt(xi^2+domeg^2);
vmat = zeros(4,4);
vmat(1,:)= [exp(-zk*y0),-exp(-zk*y0),-exp(zk*y0),0];
vmat(2,:)= [-zk*exp(-zk*y0),zk*exp(-zk*y0),-zk*exp(-zk*y0),0];
vmat(3,:) = [0,1,1,-1];
vmat(4,:) = [0,-zk,zk,-zk+2*dm];

vrhs = zeros([4,1]);
vrhs(1) = 0;
vrhs(2) = 1/sqrt(2*pi);
vrhs(3) = 0;
vrhs(4) = 0;

vcoef = vmat\vrhs;
vals(i) = vcoef(4);

end
end

rr=0.25;
x = 0.1;
y = 0;
y0 = 1;
dl = dm;

xs = tt(:).';
vs = zeros(size(xs));
for i=1:numel(xs)
x = xs(i);

f1 = @(xi) get_gxi(xi,dm,de,dl,y0,x0,x,y);
i1 = integral(f1,-15,-de-rr,'ArrayValued',true);

f2 = @(xi) (rr*1i*exp(-1i*xi))*get_gxi(-de-rr*exp(-1i*xi),dm,de,dl,y0,x0,x,y);
i2 = integral(f2,0,-pi,'ArrayValued',true);

i3 = integral(f1,-de+rr,de-rr,'ArrayValued',true);

f2 = @(xi) (rr*1i*exp(-1i*xi))*get_gxi(de-rr*exp(-1i*xi),dm,de,dl,y0,x0,x,y);
i4 = integral(f2,0,pi,'ArrayValued',true);

i5 = integral(f1,de+rr,15,'ArrayValued',true);

%i2 = 0;
%i4 = 0;
vs(i) = i1+i2+i3+i4+i5;
end

dl  = 0;
vs0 = zeros(size(xs));

for i=1:numel(xs)
x = xs(i);

f1 = @(xi) get_gxi(xi,dm,de,dl,y0,x0,x,y);
i1 = integral(f1,-15,-de-rr,'ArrayValued',true,'RelTol',10^(-14));

f2 = @(xi) (rr*1i*exp(-1i*xi))*get_gxi(-de-rr*exp(-1i*xi),dm,de,dl,y0,x0,x,y);
i2 = integral(f2,0,-pi,'ArrayValued',true);

i3 = integral(f1,-de+rr,de-rr,'ArrayValued',true);

f2 = @(xi) (rr*1i*exp(-1i*xi))*get_gxi(de-rr*exp(-1i*xi),dm,de,dl,y0,x0,x,y);
i4 = integral(f2,0,pi,'ArrayValued',true);

i5 = integral(f1,de+rr,15,'ArrayValued',true,'RelTol',10^(-14));

%i2 = 0;
%i4 = 0;
vs0(i) = i1+i2+i3+i4+i5;
end