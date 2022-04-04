% 
% clc
% clear all
% close all

% V = sqrt(3)*7200;
% I = 9e6/(V*sqrt(3));
% Id = I*sqrt(2) = 589.2557
% 100 = 17% pertubation
signal_length = 200; % length of signal in seconds
dt = 1e-4; % time increment in seconds
fs = 10;
time = 0:dt:signal_length; % time samples in seconds
amplitude = 0.2; % amplitude of wave
f0 = 0.1; % Frequency of wave at t0
f1 = 0.5; % Frequency of wave at time(end)
% Define the rectangsle signal timeseries
chirpsignal = chirp(time,f0,time(end),f1)*amplitude;
chirpsignal = chirpsignal';
squarechirpsignal = sign(chirpsignal)*amplitude;
% squarechirpsignal = squarechirpsignal';
time = time';

outvar = [time, chirpsignal, squarechirpsignal];
save final.mat outvar 
dlmwrite('chirp_sig.csv',outvar)
% save('final.mat','time','chirpsignal','squarechirpsignal','-append')

a = readtable('chirp_sig.csv');
t= a{:,1};
m = a{:,2};
m1 = a{:,3};
plot(t,m)
data_input = [t,m];
data_input1 = [t,m1];
% plot(data_input(:,1),data_input(:,2))
% % 
%  spectrogram(chirpsignal,128,120,128,1/dt)
% view(-45,65); 
% colormap bone
% f = linspace(f0,f1,6001)
% plot(time,f,'linewidth', 2)
% xlabel('Time [s]')
% ylabel('Hz')
% grid on
% set(gca, 'Fontsize', 14)
% % Plot the timeseries
% figure(1)
% plot(time,squarechirpsignal,'r-',time,chirpsignal,'b-', 'linewidth', 1)
% % plot(time,chirpsignal, 'linewidth', 1)
% xlabel('Time [s]')
% ylabel('Amplitude [p.u.]')
% title('Chirp & Square Chirp Signal')
% % title('Chirp Signal')
% grid on
% set(gca, 'Fontsize', 14);
% % xlim([0 5])
% set(gcf, 'units', 'inches', 'position', [5 0.2 5 4]);
% 
% figure(2);
% subplot(2,1,1)
% obw(squarechirpsignal,1/dt)
% title('Sq. Wave Chirp Signal: 99% Occupied Bandwidth : 3.014 Hz')
% grid on;
% subplot(2,1,2)
% obw(chirpsignal,1/dt)
% title('Chirp Signal: 99% Occupied Bandwidth : 102.338 mHz')
% 
% Y = fft(squarechirpsignal);    % Disp FFT
% % Y = Y.*conj(Y)/200;
% L = length(squarechirpsignal); % t = (0:L-1)*ts; 
% P2 = abs(Y/L); P1 = P2(1:floor(L/2+1)); 
% P1(2:end-1) = 2*P1(2:end-1); freqs = fs*(0:(L/2))/L;
% figure(3)
% subplot(2,1,1)
% plot(freqs, P1,'color', 'red', 'linewidth', 1)
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [p.u.]')
% title('Freq. Spec. of Square Chirp Signal')
% grid on
% set(gca, 'Fontsize', 14);
% % xlim([0 5])
% set(gcf, 'units', 'inches', 'position', [5 0.2 5 4]);
% 
% 
% Y = fft(chirpsignal);    % Disp FFT
% % Y = Y.*conj(Y)/200;
% L = length(chirpsignal); % t = (0:L-1)*ts; 
% P2 = abs(Y/L); P1 = P2(1:floor(L/2+1)); 
% P1(2:end-1) = 2*P1(2:end-1); freqs = fs*(0:(L/2))/L;
% subplot(2,1,2)
% plot(freqs, P1,'color', 'red', 'linewidth', 1)
% xlabel('Frequency [Hz]')
% ylabel('Magnitude [p.u.]')
% title('Freq. Spec. of Chirp Signal')
% grid on
% set(gca, 'Fontsize', 14);
% % xlim([0 5])
% set(gcf, 'units', 'inches', 'position', [5 0.2 5 4]);
% 
