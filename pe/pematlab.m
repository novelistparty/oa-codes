% Downslope propagation - Line Source
% Based on Jensen and Tindle, 1987: Mode Propagation in a Wedge.
% Equations and numerical formulation follow Comp. Ocean Acoustics, ch. 6.
tic
% Water depth at the source is 100m. 
% Range is from 0 (where source is at 100m depth) to 10km.
% The wedge is of 10 degrees.
% The wedge boundary may be written as a line of z = 0.1763 * r + 100;

% The water has parameters:
    c1 =    1500; % m/s
    rho1 =  1000; % kg/m^3
    
% The bottom has parameters:
    c2 =    1700; % m/s
    rho2 =  1500; % kg/m^3
    alpha = 0.5;    % dB/wavelength
    
% Computational absorption layer:
    alphaAL = 0.01;

% Source parameters
    zs = 180;
    f0 = 25; 
    k0 = 2 * pi * f0 / c1;
    k2 = 2 * pi * f0 / c2;
    lambda0 = c1 / f0; 

% Define zmax of computational domain Eq. (6.214)
 % H is max. depth of physical domain on propg. path 
    H = 600; 
    zmax = 4 * H / 3;
    rmax = 15000;
    
% Grid Size (6.216) and (6.219)
    dz = lambda0 / 12;
    dr = dz;
    
%% Density reduced Squared index of refraction. Eq. (6.152).
r = 0.01: dr : rmax;
z = 0: dz : zmax;

[R, Z] = meshgrid(r,z);

%% Set depth of interface
nr = length(r);
d0 = 0 * r + 200;
slope.start = round(nr/3);
slope.end = round(nr*12.5/15);
d0(slope.start:slope.end) = -.027 * r(slope.start:slope.end) + 335; % like front of COA
d0(slope.end:end) = 0;
figure; plot(r/1e3,d0)

%% Grid it
[D0, ~] = meshgrid(d0,z);

% %D0 =  0.1763 * R + 100;  % Depth to interface.
% %  D0 = - 0.1763 * R + 400;  % Depth to interface.
D = (zmax - H) / 3;       % Thickness of absorption layer (for eq. 6.146)
nsqr = ones(length(z),length(r));   % NOT density-reduced yet! (for eq. 6.158)

%% Make nsqr. Third 'if' accounts for absorption layer below phys. domain
for jj = 1:length(r)
    for ii = 1:length(z)
        if z(ii) <= d0(jj) % then we are in water
            nsqr(ii,jj) = 1;
            if z(ii) < 0.1 && d0(jj) < 0.1 % no more ocean
                nsqr(ii,jj) = (c1/c2)^2 * (1 + 1j * alpha / 27.29);
            end
        elseif z(ii) > d0(jj) && z(ii) <= H
            nsqr(ii,jj) = (c1/c2)^2 * (1 + 1j * alpha / 27.29);
        elseif z(ii) > H
            nsqr(ii,jj) = (c1/c2)^2 + 1i * alphaAL * exp(-((z(ii) - zmax)/D)^2);
        else
            error('The nsqr loop should not get here...')
        end
    end
end
%% Make dr_nsqr (Density-Reduced)

% Density is given by equation (6.153)
% By MATLAB vectorization, use wedge interface depth expr. in (6.153):
L = 2/k0;
rho = 0.5*(rho2 + rho1) + 0.5*(rho2 - rho1) * tanh((Z - D0)/L);

% First and second derivatives of rho, for use in the next expressions
drho1 = (1/(2*L)) * (rho2 - rho1) .* sech((Z - D0)/L).^2;
drho2 = (-1/(L^2)) * (rho2 - rho1) .* tanh((Z - D0)/L) .* sech((Z - D0)/L).^2;

% Compute density-reduced squared index of refraction
dr_nsqr = nsqr + (1 / (2*k0^2)) * ((drho2./rho) - 3 * drho1.^2 ./ (2*rho.^2));

% Plot of Density of field
figure; pcolor(r,z,rho)
shading interp
grid off
colorbar
set(gca,'Ydir','reverse')
% axis equal
axis tight

%% Plot of Density Reduced n^2 
% figure; 
% subplot(2,1,1); pcolor(r,z,real(dr_nsqr))
% shading interp; grid off; colorbar; set(gca,'Ydir','reverse'); axis tight
% subplot(2,1,2); pcolor(r,z,imag(dr_nsqr))
% shading interp; grid off; colorbar; set(gca,'Ydir','reverse'); axis tight

%% Starting fields - psi(r,z)
psi = ones(length(z), length(r));

%  % Greene's Source - eq. (6.102)
% psi = zeros(length(z), length(r));
% psi(:,1) = sqrt(k0) * (1.4467 - 0.4201 .* k0^2 .* (z - zs).^2)...
%         .* exp(-((k0^2 * (z - zs).^2) / 3.0512));

% Standard Gaussian Source
% psi(:,1) = sqrt(k0) * exp(-0.5 * k0^2 * (Z(:,1) - zs).^2) ... 
%             - sqrt(k0) * exp(-0.5 * k0^2 * (Z(:,1) + zs).^2);
% psi(:,1) = sqrt(k0) * exp(-0.5 * k0^2 * (Z(:,1) - zs).^2) ;

% Generalized Gaussian Source
tht1 = pi/4;
tht2 = pi/6;
psi(:,1) = sqrt(k0) * tan(tht1) * exp(-0.5*k0^2*(Z(:,1) - zs).^2*tan(tht1)^2) ...
    .* exp(1i * k0 * (Z(:,1) - zs) * sin(tht2)) ...
    - sqrt(k0) * tan(tht1) * exp(-0.5*k0^2*(Z(:,1) + zs).^2*tan(tht1)^2) ...
    .* exp(1i*k0*(Z(:,1) - 220)*sin(tht2));

% % Line Source apprx. with depth dep. of mode 3
% ind = find(z>100);
% psi(:,1) = sin(3* pi/100 * z);
% psi(:,1) = abs(psi(:,1)) + 1i * sign(psi(:,1));
% psi(ind:end,1) = 0;

% % Line Source apprx. with depth dep. of mode 2
% H0 = 410;
%  psi(:,1) = cos(2*pi/(H0) * z.' );
% psi(:,1) = abs(psi(:,1)) ;
% ind = find(z>H0);
% psi(ind:end,1) = 0;

% Make the initial field density-reduced.
dr_psi = psi;
dr_psi = dr_psi ./ sqrt(rho);       
dr_psi(1,:) = 0;

% % Plot source initial field
% figure;
% subplot(2,1,1); plot(z, real(psi(:,1)));
% subplot(2,1,2); plot(z, imag(psi(:,1)));
% % xlim([0, 150]); 
% title('Source field')
% legend('dr_\psi', 'rho')

%% Split-step Algorithm

% First apply eq. (6.129)
% N = 2^nextpow2(length(z));
N = length(z);
kzmax = 2*pi/dz;
kz = linspace(0,kzmax,N);
kz = kz.';

step = 0 * dr_psi(:,1);    
stepi = step;   % To be used for real part 
stepr = stepi;  % To be used for imag part

nR = length(r) - 1;
for jj = 1 : (length(r) - 1)
    step = exp(0.5 * 1i * k0 .* (dr_nsqr(:,jj) - 1) * dr) .* dr_psi(:,jj);
    
    stepr = dst(real(step));
    stepi = dst(imag(step));
           
    step = exp(- 1i * dr * kz.^2 / (2 * k0)) .* (stepr + 1i * stepi);
    
    stepr = idst(real(step));
    stepi = idst(imag(step));
    
    dr_psi(:,jj+1) =  exp(0.5 * 1i * k0 .* (dr_nsqr(:,jj) - 1) * dr) .* (stepr + 1i * stepi);

    disp(nR - jj)
end

% Convert density-reduced pressure to pressure 
psi = dr_psi .* sqrt(rho);
% psi = dr_psi;
% Remove absorption layer for plotting
ind = find(z>H, 1);
psi = psi(1:ind, :);
R = R(1:ind, :);
Z = Z(1:ind, :);

% %% Save Data to txt Files for Plotting in PyLab
% psiSave = -20 * log10( abs(psi) ./ sqrt(R));
% 
% save('upslope\p.txt', 'psiSave', '-ascii', '-double')
% save('upslope\R.txt', 'R', '-ascii', '-double')
% save('upslope\Z.txt', 'Z', '-ascii', '-double')
% 
% cd 'upslope'
% system('python U:\classes\ECE576_compmethods\pe_model\upslope\contour_plotting.py')
% cd ..


%% Plot TL
figure; imagesc(r/1e3, z, -20 * log10( abs(psi) ./ sqrt(R) ))
shading interp
grid off
colorbar
% axis equal
% axis tight
% xlim([0 500])
ylim([0 600])
set(gca,'Ydir','reverse')
caxis([30 100])
colormap(flip(jet))

hold on
plot(r/1e3, d0,'LineWidth',2)
plot(r(1)/1e3, zs, 'kx','LineWidth',2,'MarkerSize',10)
hold off

xlabel('Range (km)')
ylabel('Depth (m)')
title(['TL (dB), ',num2str(f0),' Hz, ',num2str(zs),' m Source'])


%% Plot TL Contours
 v = [18 22 26 30 34 38 42 46];
figure; 
contourf(R/1e3, Z, -20 * log10( abs(psi) ./ sqrt(R) ),v,'CDataMapping','direct')
colorbar
% axis equal
% axis tight
% xlim([0 500])
 ylim([0 600])
colormap('spring')
set(gca,'Ydir','reverse')

%%
toc