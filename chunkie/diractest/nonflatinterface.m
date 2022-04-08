function [r,d,d2] = nonflatinterface(t)
   xs = t;
   xp = ones(size(t));
   xpp = zeros(size(t));
   
   
    ys = exp(-t.^2/2).*sin(t);
    yp = exp(-t.^2/2).*(cos(t) - t.*sin(t));
    ypp = exp(-t.^2/2).*(-2*t.*cos(t) + (t.^2-2).*sin(t));
    
    r = [(xs(:)).'; (ys(:)).'];
    d = [(xp(:)).'; (yp(:)).'];
    d2 = [(xpp(:)).'; (ypp(:)).'];
end