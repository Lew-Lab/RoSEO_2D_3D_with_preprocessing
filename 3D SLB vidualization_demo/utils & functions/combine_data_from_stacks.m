
%% 
R_save=[];

indx = abs(SM_est_save_all(:,2))<1400  & abs(SM_est_save_all(:,3))<1400 & SM_est_save_all(:,4)>200 & SM_est_save_all(:,4)<900  & SM_est_save_all(:,5)>800;
%indx(:)=1;
r=30;
x1 = SM_est_save_all(indx,2);
y1 = SM_est_save_all(indx,3);
z1 = SM_est_save_all(indx,4);
[R,~]=sphereFit([x1,y1,z1])
%R = [0,0,0];
%
%indx = abs(SM_est_save_all(:,2))<1400  & abs(SM_est_save_all(:,3))<1400 & SM_est_save_all(:,4)>200 & SM_est_save_all(:,4)<900  & SM_est_save_all(:,5)>800;
%indx(:)=1;
r=30;
x = SM_est_save_all(:,2)-R(1);
y = SM_est_save_all(:,3)-R(2);
z = SM_est_save_all(:,4);
thetaD = Angle_save(:,2);
phiD = Angle_save(:,3);
omega = Angle_save(:,6);
Angle_save_final = Angle_save;
SM_est_final = SM_est_save_all;
%
R_save = [R_save;R];
%save('est_retrieval_1.2_v2.mat','x','y','z','thetaD','phiD','omega','Angle_save_final','SM_est_final');
%%

indx = abs(SM_est_save_all(:,2))<1400  & abs(SM_est_save_all(:,3))<1400 & SM_est_save_all(:,4)>200 & SM_est_save_all(:,4)<900  & SM_est_save_all(:,5)>800;
%indx(:)=1;
r=30;
x1 = SM_est_save_all(indx,2);
y1 = SM_est_save_all(indx,3);
z1 = SM_est_save_all(indx,4);
[R,~]=sphereFit([x1,y1,z1])
%R = [0,0,0];

r=30;
indx = abs(SM_est_save_all(:,2))<1300  & abs(SM_est_save_all(:,3))<1300 & SM_est_save_all(:,4)>100 & SM_est_save_all(:,4)<900  & SM_est_save_all(:,5)>800;
indx(:)=1;
x = [x;SM_est_save_all(:,2)-R(1)];
y = [y;SM_est_save_all(:,3)-R(2)];
z = [z;SM_est_save_all(:,4)];
thetaD = [thetaD;Angle_save(:,2)];
phiD = [phiD;Angle_save(:,3)];
omega = [omega;Angle_save(:,6)];
Angle_save_final = [Angle_save_final;Angle_save];
SM_est_final = [SM_est_final;SM_est_save_all];

R_save = [R_save;R];