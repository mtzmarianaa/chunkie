function [r,d,d2] = flatinterface(t);
   xs = t;
   xp = ones(size(t));
   xpp = zeros(size(t));
   
   
    ys = zeros(size(t));
    yp = ys;
    ypp = ys;
    
    r = [(xs(:)).'; (ys(:)).'];
    d = [(xp(:)).'; (yp(:)).'];
    d2 = [(xpp(:)).'; (ypp(:)).'];
end