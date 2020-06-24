function [taugrid,BS_start,WLS_start,qreg_start,qreg,BS_start_smooth,WLS_start_smooth] = loadsmoothbootstrap(year)

betahatsmoothBS = csvread(sprintf('BS_WLS%dsmooth.csv',year),1,0);

taugrid = betahatsmoothBS(:,1);
BS_start = betahatsmoothBS(:,2);
WLS_start = betahatsmoothBS(:,3);
qreg_start = betahatsmoothBS(:,4);
qreg = betahatsmoothBS(:,5);
BS_start_smooth = betahatsmoothBS(:,6);
WLS_start_smooth = betahatsmoothBS(:,7);
