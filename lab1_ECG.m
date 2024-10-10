%% Ali Khosravipour 99101502 - Mohamad Hosein Faramarzi 99104095
clc;
clear;
load("ECG_sig.mat");

%% Part 1
trans_Sig = Sig';
real_t = TIME(1, 1:size(trans_Sig,2));
figure;
plot(real_t, trans_Sig(1, :))
hold on;
plot(real_t, trans_Sig(2, :))
xlabel("Time(s)");
ylabel("Voltage");
title("Channels over time");
legend({'Lead I','Lead II'});

%%
figure;
numSamples = 1.5 * sfreq; 
plot(real_t(1, 1:numSamples), trans_Sig(1, 1:numSamples))
hold on;
plot(real_t(1, 1:numSamples), trans_Sig(2, 1:numSamples))
grid on;
xlabel("Time(s)");
legend({'Lead I','Lead II'});

%% Part 2

arrhythmiaLabels = {
    'NOTQRS', 'NORMAL', 'LBBB', 'RBBB', 'ABERR', 'PVC', 'FUSION', ...
    'NPC', 'APC', 'SVPB', 'VESC', 'NESC', 'PACE', 'UNKNOWN', 'NOISE', ...
    'ARFCT','','','STCH', 'TCH', 'SYSTOLE', 'DIASTOLE', 'NOTE', 'MEASURE', ...
    'PWAVE', 'BBB', 'PACESP', 'TWAVE', 'RHYTHM', 'UWAVE', 'LEARN', ...
    'FLWAV', 'VFON', 'VFOFF', 'AESC', 'SVESC', 'LINK', 'NAPC', 'PFUS', ...
    'WFON', 'WFOFF', 'RONT'
};

figure;
plot(real_t, trans_Sig(1, :)); 
hold on;
plot(real_t, trans_Sig(2, :)); 

for i = 1:length(ATRTIMED)
    rPeak = ATRTIMED(i);
    arrhythmiaType = ANNOTD(i);
    label = arrhythmiaLabels{arrhythmiaType+1};
    idx = find(real_t == rPeak);
    text(idx, trans_Sig(1, idx), label);
end


%% Part 4
fs = 360; 

% First section of ECG signal
start_time_normal = find(real_t == ATRTIMED(2));
end_time_normal = find(real_t == ATRTIMED(6));
shift = 15;
signal1 = trans_Sig(1, start_time_normal+shift:end_time_normal-shift);  % First signal
signal2 = trans_Sig(2, start_time_normal+shift:end_time_normal-shift);  % Second signal

% Plot time domain signals
figure;
plot(start_time_normal+shift:end_time_normal-shift, signal1);
hold on;
plot(start_time_normal+shift:end_time_normal-shift, signal2);
grid on;
title('Time-Domain Signals for First Section');
xlabel('Samples');
ylabel('Amplitude');

% FFT of the first signal
n1 = length(signal1);  
f1 = (0:n1-1)*(fs/n1);  
fft_signal1 = abs(fft(signal1));

figure;
plot(f1, fft_signal1);
title('FFT of First Signal (First Section)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Spectrogram of the first signal
figure;
spectrogram(signal1, 128, 120, 128, fs, 'yaxis');
title('Spectrogram of First Signal (First Section)');
xlabel('Time');
ylabel('Frequency (Hz)');


% Second section of ECG signal
start_time_normal = find(real_t == ATRTIMED(555));
end_time_normal = find(real_t == ATRTIMED(558));
signal1 = trans_Sig(1, start_time_normal+shift:end_time_normal-shift);  % First signal
signal2 = trans_Sig(2, start_time_normal+shift:end_time_normal-shift);  % Second signal

% Plot time domain signals
figure;
plot(start_time_normal+shift:end_time_normal-shift, signal1);
hold on;
plot(start_time_normal+shift:end_time_normal-shift, signal2);
grid on;
title('Time-Domain Signals for Second Section');
xlabel('Samples');
ylabel('Amplitude');

% FFT of the second signal
n2 = length(signal1);  
f2 = (0:n2-1)*(fs/n2);  
fft_signal2 = abs(fft(signal1));  

figure;
plot(f2, fft_signal2);
title('FFT of First Signal (Second Section)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Spectrogram of the second signal
figure;
spectrogram(signal1, 128, 120, 128, fs, 'yaxis');
title('Spectrogram of First Signal (Second Section)');
xlabel('Time');
ylabel('Frequency (Hz)');
