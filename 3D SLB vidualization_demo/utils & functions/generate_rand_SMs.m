function [x,y,z,thetaD,phiD] = generate_rand_SMs(n_SMs)

% generate random angular combination from uniformly sampled space
n_SMs_large = n_SMs;
x1 = rand(n_SMs_large,1)*2-1;
x2 = rand(n_SMs_large,1)*2-1;
x = 2*x1.*sqrt(1-x1.^2-x2.^2);
y = 2*x2.*sqrt(1-x1.^2-x2.^2);
z = 1-2*(x1.^2+x2.^2);

indx = z>0 | x1.^2+x2.^2>1;
x(indx)=[];
y(indx)=[];
z(indx)=[];

thetaD = atan2(sqrt(x.^2+y.^2),z)/pi*180;
phiD = atan2(y,x)/pi*180+180;

phiD(thetaD>90)=phiD(thetaD>90)+180; phiD = rem(phiD,360);
thetaD(thetaD>90)=180-thetaD(thetaD>90);
%phiD =phiD+180;
phiD(phiD>180)=-360+phiD(phiD>180);


end