function [X_ref,Y_ref,Z_ref,thetaD1,phiD1]=sphereLine(x_cord,y_cord,z_cord,R)


[X,Y,Z]=sphere(150);
X=X(:)*R; Y = Y(:)*R; Z = Z(:)*R+R;

distance = (x_cord.'-X).^2+(y_cord.'-Y).^2+(z_cord.'-Z).^2;
[~,indx] = min(distance,[],1);
X_ref = X(indx);
Y_ref = Y(indx);
Z_ref = Z(indx);

thetaD1 = real(acos((1000-(Z_ref))./1000)/pi*180); phiD1 = atan2(Y_ref,X_ref)/pi*180; 
phiD1(thetaD1>90)=phiD1(thetaD1>90)+180; thetaD1(thetaD1>90)=180-thetaD1(thetaD1>90);
phiD1 = rem(phiD1,360); phiD1(phiD1>180)=-360+phiD1(phiD1>180);


end