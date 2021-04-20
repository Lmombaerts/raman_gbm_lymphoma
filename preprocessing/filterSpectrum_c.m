function noiseFiltered_spectrum = filterSpectrum_c(raman_patientsInfo, background) 
% Filtering noise from the Raman processed counts from outliers data. 
% Outliers can originate from the "Intrinsic", "Instrument Response", 
% "background_raman_scan" or any "RawIndividualX". 
% @LaurentMombaerts 20/01/2021

% Note from manufacturer : 
% "It is assumed that the signal is further choped, because there were no
% boundary conditions applied to the pre-filtering". It also looks as if they
% have preprocessed the raman counts for the fluorescence and the
% instrument response but not the cosmic rays. The crop needs to be done at
% half the size of the SG filter + 1 to avoid boundary effects then.

%% Raman Signal

% Background scan is now the median of all backgrounds
% Instrument response is taken from data (always the same)
instrumentResp = raman_patientsInfo{1,2}.InstrumentResponse;

% RawIndividual Spectrums
rawIndividual_columnMask = startsWith(raman_patientsInfo{1,2}.Properties.VariableNames,'RawIndividual'); 
rawIndividuals = raman_patientsInfo{1,2}(:,rawIndividual_columnMask); 

% Raman shift frequencies (x-axis)
ramanshift = raman_patientsInfo{1,2}.Ramanshift1cm;

% Parameters
parameters = raman_patientsInfo{1,1}.PlaneInstanceProperties{end};

% Cropping
upper_wavelength_crop = parameters.upper_wavelength_crop;
lower_wavelength_crop = parameters.lower_wavelength_crop;
lower_freq_idx_crop = min(find(raman_patientsInfo{1,2}.Ramanshift1cm > lower_wavelength_crop));
upper_freq_idx_crop = max(find(raman_patientsInfo{1,2}.Ramanshift1cm < upper_wavelength_crop));
freqs_bin = lower_freq_idx_crop:upper_freq_idx_crop;

%% Cleaning of the Raman spectrum

% Removal of outliers. This is hand-made function that can be changed. Not
% sure whether "intrinsic_clean" is useful.
rawIndividuals_clean = removeOutliers_b(rawIndividuals);

% This is the formula reversed from the manufacturer C#/C++ code
input_filter = mean(instrumentResp .* (rawIndividuals_clean{:,:} - background'),2);
noiseFiltered_spectrum = openSG(input_filter',ramanshift,freqs_bin);

figure; subplot(211); hold on; plot(rawIndividuals_clean{:,:}); title(raman_patientsInfo{1,1}.PatientId); 
subplot(212); plot(noiseFiltered_spectrum); title(raman_patientsInfo{1,1}.ScanId); 
end