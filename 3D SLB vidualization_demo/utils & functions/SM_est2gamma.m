function [gamma,loc] = SM_est2gamma(SM_est,imgPara)

pix_sizex = imgPara.pix_sizex;
pix_sizez = imgPara.pix_sizez;

x = SM_est(:,2);
y = SM_est(:,3);
z = SM_est(:,4);
s = SM_est(:,5);
sXX = s.*SM_est(:,6);
sYY = s.*SM_est(:,7);
sZZ = s.*SM_est(:,8);
sXY = s.*SM_est(:,9);
sXZ = s.*SM_est(:,10);
sYZ = s.*SM_est(:,11);

loc = [round(x/pix_sizex)*pix_sizex,round(y/pix_sizex)*pix_sizex,round(z/pix_sizez)*pix_sizez].';
dx = (x-loc(1,:).')/100;
dy = (y-loc(2,:).')/100;
dz = (z-loc(3,:).')/100;

sXXdx = sXX.*dx;
sYYdx = sYY.*dx;
sZZdx = sZZ.*dx;
sXXdy = sXX.*dy;
sYYdy = sYY.*dy;
sZZdy = sZZ.*dy;
sXXdz = sXX.*dz;
sYYdz = sYY.*dz;
sZZdz = sZZ.*dz;

gamma = [sXX.';sYY.';sZZ.';sXY.';sXZ.';sYZ.';sXXdx.';sYYdx.';sZZdx.';sXXdy.';sYYdy.';sZZdy.';sXXdz.';sYYdz.';sZZdz.'];
end