function detrended_spectrum = openSG(meanRamanSpectrum,ramanshift_init,freqs_bin)
% Removes the polynomial trend out of the mean spectrum
% @LaurentMombaerts 04/02/2021

% There are 2 regions in the spectrum (finger, ch)
fingerWindow = 81; fingerDegree = 2;
chFullWindow = 85; chDegree = 1;
initialWindow = 11; initialDegree = 2;
openingWindow = 65; % Window size for opening operator
breakpoint = 1200; % The channel at which we transition from finger to ch.
maxIterations = 450; % Max number of iterations

% Initial Smooth signal
smoothed = smoothdata(meanRamanSpectrum,'sgolay',initialWindow,'degree',initialDegree,'samplepoints',ramanshift_init);

% Erosion of spectrum (local minimum operator)
halfWindow = (openingWindow - 1)/2;
erodedRaman = smoothed;
for i = 1:length(meanRamanSpectrum)
    % No boundary conditions since signal will be cropped 
   if (i > halfWindow) && (i < (length(meanRamanSpectrum) - halfWindow)) 
       erodedRaman(i) = min(smoothed(i-halfWindow:i+halfWindow));
   elseif i <= halfWindow
       erodedRaman(i) = min(smoothed(i:i+halfWindow));
   elseif i >= length(meanRamanSpectrum) - halfWindow
       erodedRaman(i) = min(smoothed(i-halfWindow:i));
   end
end
    
openedRaman = smoothed;
for i = 1:length(meanRamanSpectrum)
    % No boundary conditions since signal will be cropped 
   if (i > halfWindow) && (i < (length(meanRamanSpectrum) - halfWindow)) 
       openedRaman(i) = max(erodedRaman(i-halfWindow:i+halfWindow));
   elseif i <= halfWindow
       temp = max(erodedRaman(i:i+halfWindow));
       if erodedRaman(i) > temp
          openedRaman(i) = temp;
       end
   elseif i >= length(meanRamanSpectrum) - halfWindow
       temp = max(erodedRaman(i-halfWindow:i));
       if erodedRaman(i) > temp
          openedRaman(i) = temp;
       end
   end
end

% Set initial guess to the opening operator acting on the smoothed signal
tempSignal = openedRaman;
tempFinger = tempSignal(1:breakpoint);
tempCH = tempSignal(breakpoint+1:end);

localBaseline = zeros(1,length(tempSignal));
for i = 1:maxIterations
    localBaseline(1:breakpoint) = smoothdata(tempFinger,'sgolay',fingerWindow,'degree',fingerDegree,'samplepoints',ramanshift_init(1:breakpoint));
    localBaseline(breakpoint+1:end) = smoothdata(tempCH,'sgolay',chFullWindow,'degree',chDegree,'samplepoints',ramanshift_init(breakpoint+1:end));
    tempSignal = min(localBaseline,smoothed);
    
    tempFinger = tempSignal(1:breakpoint);
    tempCH = tempSignal(breakpoint+1:end);
end

% Output
detrended_spectrum = meanRamanSpectrum(freqs_bin) - localBaseline(freqs_bin);

% Plots double check
figure; x1 = subplot(211); hold on; grid; box on; 
title('Savitzky golay + sequential morphological opening fitting');
plot(ramanshift_init, meanRamanSpectrum); 
plot(ramanshift_init, localBaseline);
xlabel('Raman frequency'); ylabel('Raman intensity');
legend('Original spectrum','localBaseline');
x2 = subplot(212); hold on; grid; box on; 
plot(ramanshift_init(freqs_bin), detrended_spectrum);
title('Result of detrending algorithm');
xlabel('Raman frequency'); ylabel('Raman intensity');
linkaxes([x1,x2],'x');

figure; plot(ramanshift_init, meanRamanSpectrum); hold on;
plot(ramanshift_init, openedRaman); legend('Original spectrum','morphological opening'); 
end