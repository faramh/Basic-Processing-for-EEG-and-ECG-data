%% Load dataset
% Mohamad Hosein Faramarzi - Ali Khosravipour

load("EEG_sig.mat");
%% Q1

Xserie=Z(5,:);

ChannelNames=des.channelnames;
ChannelFive=ChannelNames(1,5);


samplingRate = 256; % in Hz, change based on your actual data
timeVector = (0:length(Xserie)-1) / samplingRate;

% Plotting
figure;
plot(timeVector, Xserie, 'LineWidth', 1.5); % Set line width for a smooth look

grid on;

xlabel('Time (s)', 'FontSize', 12);
ylabel('Amplitude (\muV)', 'FontSize', 12);
title(['EEG Channel: ' ChannelFive], 'FontSize', 14);

set(gca, 'FontSize', 10); 
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6); %
box off; 

axis tight;

saveas(gcf, 'EEG_Timeseries_Plot.png');

%% Q2

samplingRate = 256; % Frequency rate in Hz
timeVector = (0:length(Xserie)-1) / samplingRate;

ranges = [0 15; 18 40; 45 50; 50 timeVector(end)];

indices = round(ranges * samplingRate) + 1;

figure;

for i = 1:size(ranges, 1)
    rangeData = Xserie(indices(i,1):indices(i,2));
    rangeTime = timeVector(indices(i,1):indices(i,2));
    
    subplot(4, 1, i);
    plot(rangeTime, rangeData, 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Amplitude (\muV)', 'FontSize', 10);
    title(['EEG Channel: ' ChannelFive ' - Time ' num2str(ranges(i,1)) '-' num2str(ranges(i,2)) 's'], 'FontSize', 12);
    

    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
    box off;
    axis tight;
end


saveas(gcf, 'EEG_Timeseries_Subplots.png');

%% Q3
Xserie=Z(12,:);

ChannelNames=des.channelnames;
ChannelFive=ChannelNames(1,12);


samplingRate = 256; % in Hz
timeVector = (0:length(Xserie)-1) / samplingRate;

% Plotting
figure;
plot(timeVector, Xserie, 'LineWidth', 1.5); % Set line width for a smooth look

grid on;

% Labels and Title
xlabel('Time (s)', 'FontSize', 12);
ylabel('Amplitude (\muV)', 'FontSize', 12);
title(['EEG Channel: ' ChannelFive], 'FontSize', 14);

set(gca, 'FontSize', 10); 
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6); 
box off; 

axis tight;

saveas(gcf, 'EEG_Timeseries_Plot2.png');

ranges = [0 15; 18 40; 45 50; 50 timeVector(end)];

indices = round(ranges * samplingRate) + 1;

figure;

for i = 1:size(ranges, 1)
    rangeData = Xserie(indices(i,1):indices(i,2));
    rangeTime = timeVector(indices(i,1):indices(i,2));
    
    subplot(4, 1, i);
    plot(rangeTime, rangeData, 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Amplitude (\muV)', 'FontSize', 10);
    title(['EEG Channel: ' ChannelFive ' - Time ' num2str(ranges(i,1)) '-' num2str(ranges(i,2)) 's'], 'FontSize', 12);
    
%
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
    box off;
    axis tight;
end

%
saveas(gcf, 'EEG_Timeseries_Subplots2.png');



%% Q4

offset = max(max(abs(Z)))/3 ;
feq = 256 ;
ElecName = des.channelnames ;
disp_eeg(Z,offset,feq,ElecName) ;


%% Q5 - A
Xserie=Z(5,:);

ChannelNames=des.channelnames;
ChannelFive=ChannelNames(1,5);


% Sampling rate and time vector
samplingRate = 256; % Frequency rate in Hz
timeVector = (0:length(Xserie)-1) / samplingRate;

ranges = [2 7; 30 35; 42 47; 50 55];

indices = round(ranges * samplingRate) + 1;

figure;

for i = 1:size(ranges, 1)
    rangeData = Xserie(indices(i,1):indices(i,2));
    rangeTime = timeVector(indices(i,1):indices(i,2));
    
    subplot(4, 1, i);
    plot(rangeTime, rangeData, 'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Amplitude (\muV)', 'FontSize', 10);
    title(['EEG Channel: ' ChannelFive ' - Time ' num2str(ranges(i,1)) '-' num2str(ranges(i,2)) 's'], 'FontSize', 12);
    
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
    box off;
    axis tight;
end


saveas(gcf, 'EEG_Timeseries_Subplots.png');


%% Q5 - B

indices = round(ranges * samplingRate) + 1;
ranges = [2 7; 30 35; 42 47; 50 55];

figure;

for i = 1:size(ranges, 1)
    rangeData = Xserie(indices(i,1):indices(i,2));
    
    % Perform DFT 
    N = length(rangeData); % Number of samples in the range
    fftResult = fft(rangeData);
    
    % Frequency axis
    f = (0:N-1) * (samplingRate / N); % Create frequency axis
    
    P2 = abs(fftResult/N); % Two-sided spectrum
    P1 = P2(1:floor(N/2+1)); % Single-sided spectrum
    f = f(1:floor(N/2+1));   % Single-sided frequency axis
    
    subplot(4, 1, i);
    plot(f, P1, 'LineWidth', 1.5);
    grid on;
    xlabel('Frequency (Hz)', 'FontSize', 10);
    ylabel('Magnitude', 'FontSize', 10);
    title(['Frequency Spectrum for EEG Channel: ' ChannelFive ' - Time ' num2str(ranges(i,1)) '-' num2str(ranges(i,2)) 's'], 'FontSize', 12);
    
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
    box off;
    axis tight;
end

saveas(gcf, 'EEG_Frequency_Spectrum_Subplots.png');


%% Q7

figure;

for i = 1:size(ranges, 1)
    rangeData = Xserie(indices(i,1):indices(i,2));
    
    windowLength = 4 * samplingRate; % 4-second windows
    overlap = windowLength / 2; % 50% overlap
    nfft = 1024; % Number of FFT points for better frequency resolution
    
    % Compute the power spectral density using pwelch
    [Pxx, f] = pwelch(rangeData, windowLength, overlap, nfft, 256);
    
    subplot(4, 1, i);
    plot(f, 10*log10(Pxx), 'LineWidth', 1.5); % Convert to dB scale
    grid on;
    xlabel('Frequency (Hz)', 'FontSize', 10);
    ylabel('Power/Frequency (dB/Hz)', 'FontSize', 10);
    title(['Welch PSD for EEG Channel: ' ChannelFive ' - Time ' num2str(ranges(i,1)) '-' num2str(ranges(i,2)) 's'], 'FontSize', 12);
    
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
    box off;
    axis tight;
end

saveas(gcf, 'EEG_Pwelch_Frequency_Spectrum.png');

%% Q8

% Parameters for the spectrogram
L = 128; % Window length
Noverlap = 64; % Overlap between windows
nfft = L; % Number of FFT points
window = hamming(L); % Hamming window

figure;

for i = 1:size(ranges, 1)
    rangeData = Xserie(indices(i,1):indices(i,2));
    
    % Compute the spectrogram
    [S, F, T, P] = spectrogram(rangeData, window, Noverlap, nfft, samplingRate);
    P_smooth = imgaussfilt(abs(P), 1);  % Apply Gaussian filter for smoothing
        
    subplot(4, 1, i);
    imagesc(T, F, 10*log10(abs(P))); 
    axis xy; 
    colormap jet; 
    colorbar; 
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontSize', 10);
    title(['Spectrogram for EEG Channel: ' ChannelFive ' - Time ' num2str(ranges(i,1)) '-' num2str(ranges(i,2)) 's'], 'FontSize', 12);
    
% and axis limits
    set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.6);
    box off;
end

saveas(gcf, 'EEG_Spectrogram_Frequency_Spectrum.png');


%% Q8 test
L = 256;  
Noverlap = 192;  
nfft = L; 
window = hann(L); 

% plot the spectrogram
figure;
for i = 1:size(ranges, 1)
    rangeData = Xserie(indices(i,1):indices(i,2));
    
    [S, F, T, P] = spectrogram(rangeData, window, Noverlap, nfft, samplingRate);
    
    P_smooth = imgaussfilt(abs(P), 1);  % Apply Gaussian filter for smoothing
    
    subplot(4, 1, i);
    imagesc(T, F, 10*log10(abs(P))); 
    
    axis xy;
    colormap parula;  % Use 'parula' colormap for smoother color transitions
    colorbar;
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('Frequency (Hz)', 'FontSize', 10);
    title(['Smoothed Spectrogram for EEG Channel: ' ChannelFive ' - Time ' num2str(ranges(i,1)) '-' num2str(ranges(i,2)) 's'], 'FontSize', 12);
end

%% Q9
originalSamplingRate = 256;  % Original sampling rate in Hz
downsampleFactor = 2;  % Factor by which to downsample
newSamplingRate = originalSamplingRate / downsampleFactor;  % New sampling rate
intervalIndex = 2;
rangeData = Xserie(indices(intervalIndex,1):indices(intervalIndex,2));

% Low-pass filter 
nyquistFreq = newSamplingRate / 2;
[b, a] = butter(4, nyquistFreq / (originalSamplingRate / 2), 'low');  % 4th order Butterworth filter

filteredData = filtfilt(b, a, rangeData);

% Downsample the signal
downsampledData = downsample(filteredData, downsampleFactor);

downsampledTime = linspace(30, 35, length(downsampledData));

% Plot the original and downsampled signals
figure;
subplot(2, 1, 1);
plot(30:1/originalSamplingRate:35, rangeData);
title('Original Signal (30-35s)');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');

subplot(2, 1, 2);
plot(downsampledTime, downsampledData);
title(['Downsampled Signal (Factor = ' num2str(downsampleFactor) ')']);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');

N = length(downsampledData);
fftResult = fft(downsampledData);
P2 = abs(fftResult / N);
P1 = P2(1:floor(N/2+1));
f = (0:N-1) * (newSamplingRate / N);
f = f(1:floor(N/2+1));

figure;
plot(f, P1, 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('DFT of Downsampled Signal');

L = 128;  % Window length
Noverlap = 64;  % Overlap between windows
nfft = L;  % Number of FFT points
window = hamming(L);

figure;
[S, F, T, P] = spectrogram(downsampledData, window, Noverlap, nfft, newSamplingRate);
imagesc(T, F, 10*log10(abs(P))); 
axis xy;
colormap jet;
colorbar;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(['Spectrogram of Downsampled Signal (Sampling Rate = ' num2str(newSamplingRate) ' Hz)']);
