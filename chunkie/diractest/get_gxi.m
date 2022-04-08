function [val] = get_gxi(xi,dm,de,dl,y0,x0,x,y)
domeg = sqrt(dm*dm-de*de);
zk = sqrt(xi^2+domeg^2);
vmat = zeros(4,4);
vmat(1,:)= [exp(-zk*y0),-exp(-zk*y0),-exp(zk*y0),0];
vmat(2,:)= [-zk*exp(-zk*y0),zk*exp(-zk*y0),-zk*exp(zk*y0),0];
vmat(3,:) = [0,1,1,-1];
vmat(4,:) = [0,-zk,zk,-zk+2*dl];

vrhs = zeros([4,1]);
vrhs(1) = 0;
vrhs(2) = 1/sqrt(2*pi);
vrhs(3) = 0;
vrhs(4) = 0;

vcoef = vmat\vrhs;
val = vcoef(4)/sqrt(2*pi)*exp(1i*xi*(x-x0));
end

