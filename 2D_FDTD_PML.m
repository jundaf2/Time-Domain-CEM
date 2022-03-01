%% convention 
clear all;
close all;
clc;
%% preparation
epsilon0 = 8.8541878128E-12; % Vacuum permittivity
mu0 = 4*pi*1E-7; % Vacuum permeability
c = 1/sqrt(epsilon0*mu0); % Vacuum speed of light
eta0 = sqrt(mu0/epsilon0); % wave impedance

freq_0 = 1e8; % what frequency do we want to simulate
omega_0 = 2*pi*freq_0; % angular frequency
lambda = c/freq_0; % wavelength

N = 3000; % number of time steps for simulation
Dim = 256; % size of simulated square region in space
dx = lambda/30; % choose grid size according to page 390 of text book
dy = dx; % assume square grid
dt = 1/c/sqrt(1/dx^2+1/dy^2); % choose time step size according to page 395 of text book

% to save memory, we only store the data for the current and previous time step
Ez_x = zeros(2,Dim,Dim); 
Ez_y = zeros(2,Dim,Dim); 
Ez = zeros(2,Dim,Dim); 
Hx = zeros(2,Dim,Dim); 
Hy = zeros(2,Dim,Dim); 

%% time domain source
t = (0:N-1)*dt;N*dt/2
tau_p = 3/omega_0;%N*dt/8; % omega_0*tau_p = 3 --> tau_p = 3/omega_0
Source1 = exp(-0.5*((t-N*dt/4)/tau_p).^2);
Source2 = -(t-tau_p*5)/tau_p .* exp(-0.5*((t-tau_p*5)/tau_p).^2);
Source3 = exp(-0.5*((t-N*dt/4)/tau_p).^2) .* sin(omega_0*t); % 
Source4 = 1*(1-exp(-(t)/tau_p)).*sin(omega_0 * t); % single freq
% source position
src_x = floor(Dim/2);
src_y = floor(Dim/2);
figure(1);
plot(t,Source1,'linewidth',3,'DisplayName','Gaussian Pulse');hold on;
plot(t,Source2,'linewidth',3,'DisplayName','Neumann Pulse');hold on;
plot(t,Source3,'linewidth',3,'DisplayName','Modulated Gaussian Pulse');hold on;
plot(t,Source4,'linewidth',3,'DisplayName','Time Harmonic Sine with Taper');hold on;
legend;
xlabel('Time (s)');
ylabel('Magnitude');
% title('Gaussian');
set(gca,'fontsize',36);
set(gca,'linewidth',3);
set(gca, 'LooseInset', [0,0,0,0]);

% saveas(gca, 'pulses.png'); 
Source = Source2;

%% frequency domain
Fs = 1/dt; NFFT = 2048;
f = Fs*(0:(N/2))/NFFT;
Y = fft(Source1,NFFT);
P = abs(Y/N);
P1 = P(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
Y = fft(Source2,NFFT);
P = abs(Y/N);
P2 = P(1:N/2+1);
P2(2:end-1) = 2*P2(2:end-1);
Y = fft(Source3,NFFT);
P = abs(Y/N);
P3 = P(1:N/2+1);
P3(2:end-1) = 2*P3(2:end-1);
Y = fft(Source4,NFFT);
P = abs(Y/N);
P4 = P(1:N/2+1);
P4(2:end-1) = 2*P4(2:end-1);
figure(2)
plot(f,P1,'linewidth',3,'DisplayName','Gaussian Pulse');hold on;
plot(f,P2,'linewidth',3,'DisplayName','Neumann Pulse');hold on;
plot(f,P3,'linewidth',3,'DisplayName','Modulated Gaussian Pulse');hold on;
plot(f,P4,'linewidth',3,'DisplayName','Time Harmonic Sine with Taper');hold on;
ylim([0,1])
xlim([0,freq_0*2])
legend;
xlabel('frequency (Hz)')
ylabel('Magnitude')
set(gca,'fontsize',36);
set(gca,'linewidth',3);
set(gca, 'LooseInset', [0,0,0,0]);
% saveas(gca, 'pulseSpectrum.png'); 
%% geometry
sheet_left = floor(Dim*2/3);
sheet_right = floor(Dim*3/4);
slot_down = floor(Dim*52/100);
slot_up = floor(Dim*48/100);

dx*(slot_down-slot_up)

epsilon = epsilon0*ones(Dim,Dim); % permitivity
mu = mu0*ones(Dim,Dim); % permeability
sigma = zeros(Dim,Dim); % conductivity
sigma_x = zeros(Dim,Dim);
sigma_y = zeros(Dim,Dim);
% PML of length NPML at boundary grids
NPML = 64;
mPML = 2; % m = 2 or 3 is a good choice
R0 = 1e-12; % magnitude of reflection coefficient
sigma_max = -(mPML+1)/(2*eta0*NPML*dx)*log(R0); % equation on page of 416 textbook

for i=1:NPML
    sigma_x(i,:) = sigma_max * ((NPML+2-i)/NPML )^mPML;
end
for i=1:NPML
    sigma_y(:,i) = sigma_max * ((NPML+2-i)/NPML )^mPML;
end
for i=Dim-NPML+1:Dim
    sigma_x(i,:) = sigma_max * ((1+i-Dim+NPML)/NPML )^mPML;
end
for i=Dim-NPML+1:Dim
    sigma_y(:,i) = sigma_max * ((1+i-Dim+NPML)/NPML )^mPML;
end

figure(10);
[X,Y] = meshgrid((1:Dim)*dx,(1:Dim)*dy);
contourf(X,Y,sigma_x+sigma_y); 
colorbar;
xlabel('Y (m)')
ylabel('X (m)')
set(gca,'fontsize',36);
%set(gca, 'LooseInset', [0,0,0,0]);
%saveas(gca, 'pml_visualization.png'); 

% conductor object information
sigma_x(1:slot_up,sheet_left:sheet_right)=1e12;
sigma_x(slot_down:Dim,sheet_left:sheet_right)=1e12;
sigma_y(1:slot_up,sheet_left:sheet_right)=1e12;
sigma_y(slot_down:Dim,sheet_left:sheet_right)=1e12;

beta_x = epsilon(2:Dim-1,2:Dim-1)/dt + sigma_x(2:Dim-1,2:Dim-1)/2; 
beta_y = epsilon(2:Dim-1,2:Dim-1)/dt + sigma_y(2:Dim-1,2:Dim-1)/2;
alpha_x = epsilon(2:Dim-1,2:Dim-1)/dt - sigma_x(2:Dim-1,2:Dim-1)/2;
alpha_y = epsilon(2:Dim-1,2:Dim-1)/dt - sigma_y(2:Dim-1,2:Dim-1)/2;

% observing point
obsrv_x = floor(Dim*5/10);
obsrv_y = floor(Dim*9/10);
obsrv_data = zeros(1,N);

%% move data in matrices to GPU, if there is a GPU
try
    gpuArray(1);
    canUseGPU=true;
    gpuDevice
catch
    canUseGPU=false;
end
% canUseGPU=false;
% memory need in total: (Dim^2*(2*5+6)+N+2)*4/1e9
if canUseGPU
    Ez_x=gpuArray(single(Ez_x));
    Ez_y=gpuArray(single(Ez_y));
    Ez=gpuArray(single(Ez));
    Hx=gpuArray(single(Hx));
    Hy=gpuArray(single(Hy));

    beta_x=gpuArray(single(beta_x));
    beta_y=gpuArray(single(beta_y));
    alpha_x=gpuArray(single(alpha_x));
    alpha_y=gpuArray(single(alpha_y));
%     dx=gpuArray(single(dx));
%     dy=gpuArray(single(dy));
%     Source=gpuArray(single(Source));
    epsilon=gpuArray(single(epsilon));
    mu=gpuArray(single(mu));
    whos Ez_x
    whos Ez_y
    whos Ez
    whos Hx
    whos Hy
    whos beta_x
    whos beta_y
    whos alpha_x
    whos alpha_y
%     whos dx
%     whos dy
    whos epsilon
    whos mu
%     whos Source
end

%% 2D TM Polarized EM field FDTD calculation (page 397-398 of the textbook)
h3=figure(3);
h4=figure(4);
h5=figure(5);
tic;
for n = 2:N  % index of time
    Ez_x(2,2:Dim-1,2:Dim-1) = 1 ./ beta_x .* (alpha_x .* squeeze(Ez_x(2-1,2:Dim-1,2:Dim-1)) + 1/dx * squeeze(Hy(2-1,2:Dim-1,2:Dim-1)-Hy(2-1,1:Dim-2,2:Dim-1)));
    Ez_y(2,2:Dim-1,2:Dim-1) = 1 ./ beta_y .* (alpha_y .* squeeze(Ez_y(2-1,2:Dim-1,2:Dim-1)) - 1/dy * squeeze(Hx(2-1,2:Dim-1,2:Dim-1)-Hx(2-1,2:Dim-1,1:Dim-2)));
    Ez(2,2:Dim-1,2:Dim-1) =  Ez_x(2,2:Dim-1,2:Dim-1) + Ez_y(2,2:Dim-1,2:Dim-1); 
    Ez(2,src_x,src_y) =  Ez(2,src_x,src_y) - 1/beta_x(src_x,src_y) * Source(n-1);
    Hx(2,2:Dim-1,2:Dim-1) = 1 ./ beta_y .* (alpha_y .* squeeze(Hx(2-1,2:Dim-1,2:Dim-1)) - epsilon(2:Dim-1,2:Dim-1)./(mu(2:Dim-1,2:Dim-1)*dy) .* squeeze(Ez(2,2:Dim-1,3:Dim) - Ez(2,2:Dim-1,2:Dim-1)));
    Hy(2,2:Dim-1,2:Dim-1) = 1 ./ beta_x .* (alpha_x .* squeeze(Hy(2-1,2:Dim-1,2:Dim-1)) + epsilon(2:Dim-1,2:Dim-1)./(mu(2:Dim-1,2:Dim-1)*dx) .* squeeze(Ez(2,3:Dim,2:Dim-1) - Ez(2,2:Dim-1,2:Dim-1)));
    
    % E_\parallel = 0
    Ez(:,slot_up,sheet_left:sheet_right)=0;
    Ez(:,slot_down,sheet_left:sheet_right)=0;
    Ez(:,slot_down:end,sheet_left)=0;
    Ez(:,slot_down:end,sheet_right)=0;
    Ez(:,1:slot_up,sheet_left)=0;
    Ez(:,1:slot_up,sheet_right)=0;
    
    Ez_x(:,slot_up,sheet_left:sheet_right)=0;
    Ez_x(:,slot_down,sheet_left:sheet_right)=0;
    Ez_x(:,slot_down:end,sheet_left)=0;
    Ez_x(:,slot_down:end,sheet_right)=0;
    Ez_x(:,1:slot_up,sheet_left)=0;
    Ez_x(:,1:slot_up,sheet_right)=0;
    
    Ez_y(:,slot_up,sheet_left:sheet_right)=0;
    Ez_y(:,slot_down,sheet_left:sheet_right)=0;
    Ez_y(:,slot_down:end,sheet_left)=0;
    Ez_y(:,slot_down:end,sheet_right)=0;
    Ez_y(:,1:slot_up,sheet_left)=0;
    Ez_y(:,1:slot_up,sheet_right)=0;
    
    % H_\perp=0
    Hx(:,slot_up,sheet_left:sheet_right)=0;
    Hx(:,slot_down,sheet_left:sheet_right)=0;
    Hy(:,1:slot_up,sheet_left)=0;
    Hy(:,1:slot_up,sheet_right)=0;
    Hy(:,slot_down:end,sheet_left)=0;
    Hy(:,slot_down:end,sheet_right)=0;
    
    Ez_x(1,:,:)=Ez_x(2,:,:);
    Ez_y(1,:,:)=Ez_y(2,:,:);
    Ez(1,:,:)=Ez(2,:,:);
    Hx(1,:,:)=Hx(2,:,:);
    Hy(1,:,:)=Hy(2,:,:);
    
    
    obsrv_data(n) = Ez(2,obsrv_x,obsrv_y);
    
    if mod(n,50)==0
        range = max(max(squeeze(abs(Ez(2,:,:)))))/10;
        figure(3)
        imagesc((1:Dim)*dx,(1:Dim)*dy,abs(squeeze(Ez(2,:,:)))/range);%,'Interpolation','bilinear');
        colormap summer
        axis image
        xlabel('Y (m)')
        ylabel('X (m)')
        title(['$\mathcal{E}_z$', ' at ', 'time = ',num2str(round(n*dt,9)),' s',],'interpreter','latex');
        set(gca,'fontsize',36);
        set(gca,'linewidth',3);
        %set(gca, 'LooseInset', [0,0,0,0]);
        drawnow;
%         jFrame = get(h3,'JavaFrame');
%         set(jFrame,'Maximized',1);
        set(gcf,'Position',[0,0,1024,1024]); 
        pause(0.3);	
        
        saveas(gca, ['sim1step',num2str(n),'Ez','.png']);

        %range = max(max(squeeze(abs(Hx(2,:,:)))));
        figure(4)
        imagesc((1:Dim)*dx,(1:Dim)*dy,eta0*abs(squeeze(Hx(2,:,:)))/range);
        colormap summer
        axis image
        xlabel('Y (m)')
        ylabel('X (m)')
        title(['$\mathcal{H}_x$', ' at ', 'time = ',num2str(round(n*dt,9)),' s',],'interpreter','latex');
        set(gca,'fontsize',36);
        set(gca,'linewidth',3);
%         jFrame = get(h4,'JavaFrame');
%         set(jFrame,'Maximized',1);
        set(gcf,'Position',[0,0,1024,1024]); 
        pause(0.3);	
        saveas(gca, ['sim1step',num2str(n),'Hx','.png']);
        %set(gca, 'LooseInset', [0,0,0,0]);

        %range = max(max(squeeze(abs(Hy(2,:,:)))));
        figure(5)
        imagesc((1:Dim)*dx,(1:Dim)*dy,eta0*abs(squeeze(Hy(2,:,:)))/range);
        colormap summer
        axis image
        xlabel('Y (m)')
        ylabel('X (m)')
        title(['$\mathcal{H}_y$', ' at ', 'time = ',num2str(round(n*dt,8)),' s',],'interpreter','latex');
        set(gca,'fontsize',36);
        set(gca,'linewidth',3);
%         jFrame = get(h5,'JavaFrame');
%         set(jFrame,'Maximized',1);
        set(gcf,'Position',[0,0,1024,1024]); 
        pause(0.3);	
        saveas(gca, ['sim1step',num2str(n),'Hy','.png']);
        %set(gca, 'LooseInset', [0,0,0,0]);
    end
end
toc;
%% move the data from GPU to CPU, if there is a GPU
if canUseGPU
    Ez_x=gather(Ez_x);
    Ez_y=gather(Ez_y);
    Ez=gather(Ez);
    Hx=gather(Hx);
    Hy=gather(Hy);
    beta_x=gather(beta_x);
    beta_y=gather(beta_y);
    alpha_x=gather(alpha_x);
    alpha_y=gather(alpha_y);
    dx=gather(dx);
    dy=gather(dy);
    Ez_y=gather(Ez_y);
    epsilon=gather(epsilon);
    mu=gather(mu);
    Source=gather(Source);

end
%% further processing at observation point
figure(6)
plot(t,obsrv_data,'linewidth',3);
xlabel('Time (s)');
ylabel('Magnitude');
set(gca,'fontsize',36);
set(gca,'linewidth',3);
set(gca, 'LooseInset', [0,0,0,0]);

Fs = 1/dt; NFFT = 2048;
f = Fs*(0:(N/2))/N;
Y = fft(obsrv_data);
P = abs(Y/N);
P1 = P(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure(7)
plot(f,P1,'linewidth',3);
xlabel('frequency (Hz)')
ylabel('Magnitude')
set(gca,'fontsize',36);
set(gca,'linewidth',3);
set(gca, 'LooseInset', [0,0,0,0]);

