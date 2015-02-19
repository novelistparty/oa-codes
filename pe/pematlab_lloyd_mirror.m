% PE split-step Lloyd-mirror
% Equations and numerical formulation follow Comp. Ocean Acoustics, ch. 6.
clear all
tic

% The water has parameters:
    c1 =    1500; % m/s
    
% The bottom has parameters:
    c2 =    1800; % m/s
    alpha = 1;    % dB/wavelength
    
% Computational absorption layer:
    alphaAL = 0.01;

% Source parameters
    zs = 75;
    f0 = 100; 
    k0 = 2 * pi * f0 / c1;
    k2 = 2 * pi * f0 / c2;
    lambda0 = c1 / f0; 

% Define zmax of computational domain Eq. (6.214)
 % H is max. depth of physical domain on propg. path 
    H = 800; 
    zmax = 4 * H / 3;
    rmax = 500;
    
% Grid Size (6.216) and (6.219)
    dz = lambda0 / 6;
    dr = dz;
    
%% Density reduced Squared index of refraction. Eq. (6.152).

r = 0.01 : dr : rmax;
z = 0 : dz : zmax;
[R, Z] = meshgrid(r,z);

D0 = 1200;

D = (zmax - H) / 3;       % Thickness of absorption layer (for eq. 6.146)
nsqr = ones(length(z),length(r));   % NOT density-reduced yet! (for eq. 6.158)

%% Make nsqr. Third 'if' accounts for absorption layer below phys. domain
for jj = 1:length(r)
    for ii = 1:length(z)
        if z(ii) <= D0
            nsqr(ii,jj) = 1;
        elseif z(ii) > D0 && z(ii) <= H
            nsqr(ii,jj) = (c1/c2)^2 * (1 + 1j * alpha / 27.29);
        elseif z(ii) > H
            nsqr(ii,jj) = (c1/c2)^2 + 1i * alphaAL * exp(-((z(ii) - zmax)/D)^2);
        else
            error('The nsqr loop should not get here...')
        end
    end
end
%% Plot of n^2 
figure; 
subplot(2,1,1); pcolor(r,z,real(nsqr))
shading interp; grid off; colorbar; set(gca,'Ydir','reverse'); axis tight
subplot(2,1,2); pcolor(r,z,imag(nsqr))
shading interp; grid off; colorbar; set(gca,'Ydir','reverse'); axis tight
title('Squared Index of Refraction')

%% Starting fields - psi(r,z)
psi = ones(length(z), length(r));

% %  % Greene's Source - eq. (6.102)
% psi(:,1) = sqrt(k0) * (1.4467 - 0.4201 .* k0^2 .* (Z(:,1) - zs).^2)...
%         .* exp(-((k0^2 * (Z(:,1) - zs).^2) / 3.0512));

% Standard Gaussian Source
psi(:,1) = sqrt(k0) * exp(-0.5 * k0^2 * (Z(:,1) - zs).^2) ... 
            - sqrt(k0) * exp(-0.5 * k0^2 * (Z(:,1) + zs).^2);
% psi(:,1) = sqrt(k0) * exp(-k0^2 / 2 * (z - zs).^2) ;

% % Generalized Gaussian Source
% tht1 = pi/4;
% tht2 = pi/3;
% psi = zeros(length(z), length(r));
% psi(:,1) = sqrt(k0) * tan(tht1) * exp(-k0^2/2*(z - zs).^2*tan(tht1)^2) ...
%     .* exp(1i * k0 * (z - zs) * sin(tht2));
    
%% Plot source initial field
figure
plot(z,angle(psi(:,1)));
% xlim([0, 150]); 
title('Source field')
legend('\psi', 'rho')

%% Split-step Algorithm

% First apply eq. (6.129)
% N = 2^nextpow2(8 *(zmax/lambda0));
N = length(z);
kzmax = 2*pi/dz;
kz = linspace(0,kzmax,N);
kz = kz.';

step = 0 * psi(:,1);    
stepi = step;   % To be used for real part 
stepr = stepi;  % To be used for imag part 
nR = length(r) - 1;
for jj = 1 : nR
    step = exp(1i * k0 / 2 .* (nsqr(:,jj) - 1) * dr) .* psi(:,jj);
    
    stepr = dst(real(step));
    stepi = dst(imag(step));
           
    step = exp(- 1i * dr / (2 * k0) * kz.^2 ) .* (stepr + 1i * stepi);
    
    stepr = idst(real(step));
    stepi = idst(imag(step));
    
    psi(:,jj+1) =  exp(1i * k0 / 2 .* (nsqr(:,jj) - 1) * dr) .* (stepr + 1i * stepi);

    disp(num2str(nR - jj))
end
TL = -20 * log10( abs(psi) ./ sqrt(R.^2 + Z.^2));

%% Save Data to txt Files for Plotting in PyLab
% psiSave = -20 * log10( abs(psi) ./ sqrt(R));
% 
% save('psimat.txt', 'psiSave', '-ascii', '-double')
% save('R.txt', 'R', '-ascii', '-double')
% save('Z.txt', 'Z', '-ascii', '-double')
% psimat = load('psimat.txt');

%% Plot TL
figure; 
imagesc(r*1e-3, z, TL)
shading interp
grid off
colorbar
axis tight
ylim([0 H])
set(gca,'Ydir','reverse')
caxis([15 90])
colormap(flip(jet))

xlabel('Range (km)')
ylabel('Depth (m)')
title(['TL (dB), ',num2str(f0),' Hz, ',num2str(zs),' m Source'])

%% Plot TL Contours
v = [25, 31, 37, 43, 49, 55];
figure; contourf(R, Z, TL, v)
colorbar
xlim([0 500])
ylim([0 500])
colormap 'pink'
set(gca,'Ydir','reverse')

%% 
toc