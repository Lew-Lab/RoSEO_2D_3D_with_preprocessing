function gamma_new = gamma_update_2_loc_rec(loc_remove,gamma1,scaling_factor)


S_scdM = gamma1(1:6);
M_dx = gamma1(7:9);
M_dy = gamma1(10:12);
M_dz = gamma1(13:15);

M_dx_new = M_dx-S_scdM(1:3)*loc_remove(1)/scaling_factor(1);
M_dy_new = M_dy-S_scdM(1:3)*loc_remove(2)/scaling_factor(2);
M_dz_new = M_dz-S_scdM(1:3)*loc_remove(3)/scaling_factor(3);

gamma_new = [S_scdM;M_dx_new;M_dy_new;M_dz_new];

end