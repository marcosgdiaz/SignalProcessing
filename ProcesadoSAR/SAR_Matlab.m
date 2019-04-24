% SAR SIMULATION 

clear all,clc,close all;

fini = 14e9;        % Start frequency
ffin = 26e9;        % Stop frequency
df = 500*1e6;       % Sampling frequency (equivalent to "PRF")

Runamb = 3e8/(2*df);
disp(['Maximum unambiguity distance: ' num2str(Runamb*100) ' cm']);

f = [fini:df:ffin]; % Radar bandwidth
BW=(f(end)-f(1))/1e9;

disp(['Bandwidth = ' num2str(BW) ' GHz']);
disp(['Range (distance) resolution) = ' num2str(3e8/(2*BW*1e7)) ' cm']);
disp(' ');

lambda_fMax=3e8/max(f);

dx = 0.25; % Synthetic aperture sampling as a function of the wavelength
dx=lambda_fMax*dx;       

Lap = 0.4;              % Synthetic aperture size
x = [-Lap/2 : dx : Lap/2];
y = 0*x + 0.2;       % Position of the radar in the y-axis

% Targets position ------------------------------
xp = [-0.08 -0.05 -0.02 0 0.03 0.062];
yp = [ 0.56  0.52 0.54 0.57 0.59 0.52];

% RCS of the targets (dBm2)
rcs = [-3 -20 -10 0 -15 -6];
rcs = 10.^(rcs/10);

mean_yp=mean(yp);
dist_target=mean_yp - y(1);
disp(['Aperture size = ' num2str(Lap*100) ' cm, aperture electric size = ' num2str(Lap/lambda_fMax) ' wavelengths']);
disp(['Bandwidth = ' num2str(BW) ' GHz']);
disp(['Cross-range resolution (approx). = ' num2str(100*lambda_fMax*dist_target/Lap) ' cm']);


k = 2*pi.*f/3e8;

% Calculation of the scattered field
sigma = zeros(numel(f),numel(x));
sigma_log = sigma;

for q= 1:numel(f)
    for p=1:numel(x)
        Rp = sqrt((x(p)-xp).^2+(y(p)-yp).^2);
        sigma (q,p) = sum((rcs.*exp(-1*j*2*k(q)*Rp))./(Rp.^2));
        sigma_log(q,p) = 10*log10(abs(sigma(q,p)));
    end
       sigma_norm(q,:) = sigma_log(q,:)- max(sigma_log(q,:));
end

% FIGURE 1 - SCATTERED FIELDS
figure
plot(x,sigma_norm(1,:),'linewidth',2);
hold on
plot(x,sigma_norm((end+1)/2,:),'r','linewidth',2);
hold on
plot(x,sigma_norm(end,:),'k','linewidth',2);
xlabel('x(m)');
ylabel('Scattered field (dBV/m)');
legend('f = 18GHz','f = 22GHz','f = 26GHz');
axis([min(x) max(x) -20 0]);
grid on; box on;


% SAR imaging =============================================================

% Imaging domain
dL = dx/2;
X0 = [x(1) : dL: x(end)];
Y0 = [0.2 : dL: 0.6] + y(1);
[X1 Y1] = meshgrid(X0,Y0);
xprima = reshape (X1,numel(X1),1);
yprima = reshape (Y1,numel(Y1),1);

%reflectivity
refl = zeros(numel(xprima),1);

% Coherent combination of frequencies and aperture positions
for q = 1:numel(f)
     for p = 1:numel(x)
            Rp = sqrt((x(p)-xprima).^2+(y(p)-yprima).^2);
            refl =refl+ sigma(q,p).*exp(1*j*2*k(q)*Rp);
     end
end
reflXY = refl/max(abs(refl));
reflXY = reshape(reflXY,size(X1));
reflXY = 10*log10(abs(reflXY));

% FIGURE 2 - SAR IMAGE
figure;
hold on;

%Plots the synthetic aperture radar positions
plot(x,y,'bo','markerfacecolor','b');
grid on;

% Plots the reflectivity - SAR image -
pcolor ( X1,Y1, reflXY);
shading interp;
set(gcf,'renderer','zbuffer');
xlabel('x axis (m)')
ylabel('y axis (m)')
title(['SAR image, reflectivity, in dB, BW = ' num2str(BW) ' GHz']);
caxis([-40 0]);
grid on; colorbar;
set(gca,'dataaspectratio',[1 1 1]);

% Plots the targets with marker size proportional to RCS
rcs=10*log10(rcs);
rcs_size=[4:1:16];
rcs_scale=linspace(min(rcs),max(rcs),numel(rcs_size));
for n=1:numel(rcs)
    [kk ind]=min(abs(rcs(n)-rcs_scale));
    plot(xp(n),yp(n),'ko','markersize',rcs_size(ind),'markerfacecolor','none','linewidth',2);
end
box on;

