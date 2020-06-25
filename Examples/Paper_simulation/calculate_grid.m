function [taugrid_qreg, taugrid_midpoint, taugrid_ue] = calculate_grid(ntau)

taugrid_qreg = (1:ntau)/(ntau+1);
taugrid_midpoint = (1/2*([0:(ntau-1)]  + [1:(ntau)]))/ntau;
taugrid_ue = taugrid_midpoint;
taugrid_ue(1) = 0;
taugrid_ue(end) = 1;
