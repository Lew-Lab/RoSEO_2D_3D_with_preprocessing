function [theta, phi, alpha] = symmCone2angle(mux,muy,muz,gamma)
phi = (atan2d(muy,mux));
%muz=sqrt(1-mux^2-muy^2);
theta = (acosd(muz/1));
alpha = (acosd(sqrt(2*gamma+0.25)-0.5));
end