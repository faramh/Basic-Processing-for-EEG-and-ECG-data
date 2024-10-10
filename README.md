# Basic-Processing-for-EEG-and-ECG-data
Basic Processing for EEG and ECG data using 

Extracting the Signal:

We selected a specific time intervals from the EEG signal to analyze and manipulate further.
Low-Pass Filtering:

Before downsampling, we applied a 4th-order Butterworth low-pass filter. This filter ensures that no frequency content above half the new sampling rate (Nyquist frequency) passes through. This step prevents aliasing, a phenomenon where higher frequencies fold back into the lower spectrum during downsampling.
Downsampling:

After filtering, we downsampled the signal by a factor of 2. The original sampling rate was 256 Hz, so the new sampling rate became 128 Hz.
We created a new time vector to match the downsampled signal and plotted both the original and downsampled signals for comparison.
Discrete Fourier Transform (DFT):

We calculated the DFT of the downsampled signal to analyze the frequency content. This showed us the frequency components of the downsampled signal and allowed us to compare the magnitude of the frequency spectrum with the original signal's DFT.
Short-Time Fourier Transform (STFT / Spectrogram):

We computed the spectrogram (time-frequency representation) of the downsampled signal using the STFT with a Hamming window of length 128, a 50% overlap, and 128 FFT points.
The resulting spectrogram was plotted to visualize the time-varying frequency content of the downsampled signal, limited to the new Nyquist frequency (half of the new sampling rate).
