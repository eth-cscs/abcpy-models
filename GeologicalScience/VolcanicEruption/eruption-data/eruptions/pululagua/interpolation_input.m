wtpct = [0.139, 0.4, 1.65, 3.27, 5.96, 9.346, 12.756, 16.72, 21.037, 18.871, 7.089, 1.082, 0.603, 0.716, 0.26, 0.082, 0.016, 0.003];
phi = linspace(-7, 10, 18);
density = [600, 600, 600, 600, 600, 600, 600, 837.5, 1075,1312.5, 1550, 1787.5, 2025, 2262.5, 2500, 2500, 2500, 2500];

bar(phi, wtpct)
figure
plot(phi, density)
figure


phi_interp = linspace(-7, 10, 180);
wtpct_interp = interp1(phi, wtpct, phi_interp);
%normalize to have sum = 100
wtpct_interp = wtpct_interp ./ (sum(wtpct_interp)/100);
density_interp = interp1(phi, density, phi_interp);

bar(phi_interp, wtpct_interp)
figure
plot(phi_interp, density_interp)

pause