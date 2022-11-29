function rez = plotCalibrationInfo(varargin)
%PLOTCALIBRATIONINFO
if nargin < 1
    pxsize_um = 7.5;
else
    pxsize_um = varargin{1};
end
pstr = getPathsFromUser();

specieslist = {'mouse_g','mouse_s', 'marmoset'};
Nspecies = numel(specieslist);

fw = 35; fh = 20;
f = figure('Color','w','Units', 'centimeters');
f.ToolBar = 'none'; f.MenuBar ='none';
f.Position = [3 3 fw fh]; 

p = panel();

p.pack('v', {0.4 0.6})
p(1).pack('h', 4);
p(2).pack('h', Nspecies);

p.de.margin = 1;
p(1).de.marginleft = 20;
p(2).margintop = 25;
p(2).de.marginleft = 20;

p.margin = [15 12 6 20];

%% load the projector spectrum
dataprojspectrum      = load(pstr.dpspectrum);
dataprojspectrumblank = load(pstr.dpspectrumblank);

data_wavelengths = dataprojspectrum(:, 1);
dw_nm            = mean(diff(data_wavelengths));
projspectrum     = dataprojspectrum(:, 2) - dataprojspectrumblank(:, 2); % remove background

% setup relevant wavelength limits
start_nm     = 300; 
end_nm       = 700;
wavelengths  = start_nm:dw_nm:end_nm;

% resample projector spectrum to make sure bins are equal
projspectrum = interp1(data_wavelengths, projspectrum, wavelengths);

% denoise projector spectrum and normalize to sum
projspectrum = smoothdata(projspectrum, 'gaussian', 1/dw_nm);
projspectrum(projspectrum < 0.001) = 0;
projspectrum = projspectrum/sum(projspectrum); 


p(1, 1).select();
xlim([start_nm end_nm]);
ylim([0 max(projspectrum)])
line(wavelengths, projspectrum, 'Color', 'k')
xlabel('Wavelength (nm)')
ylabel('Power output (norm.)')
title('Device spectrum (corrected)')

%% load diode calibration 
datadiodecalibration = load(pstr.dpdiodesens);
if isstruct(datadiodecalibration)
    datadiodecalibration = datadiodecalibration.calibrationData;
end
diodsens_AperWpernm  = interp1(datadiodecalibration(:,1), datadiodecalibration(:,2), wavelengths);
p(1, 2).select();
xlim([start_nm end_nm]); ylim([0 max(diodsens_AperWpernm)])
line(wavelengths, diodsens_AperWpernm, 'Color', 'k')
title('Photodiode sensitivity')
ylabel('Sensitivity (A/W)')
xlabel('Wavelength (nm)')
%% load apperture data 
Gid = fopen(pstr.dpaperturetable);
apptabledat = textscan(Gid,'%f\t%u\t%f\t%s');
fclose(Gid);

appcurr_nA   = scalePhotodiodeToNanoAmps(apptabledat{3}, apptabledat{4});
appsizes_mm2 = (double(apptabledat{2}) * pxsize_um * 1e-3).^2;
pfit         = polyfit(appsizes_mm2, appcurr_nA, 1);

current_per_area_nApermm2 = pfit(1);

% get average irradiance (weigthed by the projector spectrum)
totalSens_AperW    = diodsens_AperWpernm * projspectrum';
irradiance_Wpermm2 = current_per_area_nApermm2 * 1e-9 / totalSens_AperW;

% transform to photons mm–2 s–1 nm–1

hplanck  = 4.135667e-15; % Planck's constant eV*s
clight   = 299792458;    % c in m/s
eV_per_J = 6.242e18;     % eV per J

% unfold average irradiance over the spectrum
irradiance_Wpermm2_per_nm = irradiance_Wpermm2 * projspectrum;

% divide power by single photon energy
photons_per_smm2_per_nm    = (irradiance_Wpermm2_per_nm * eV_per_J/(hplanck * clight)) ...
    .* (1e-9 * wavelengths);

% number of photons reaching the array
photon_flux = sum(photons_per_smm2_per_nm);

p(1, 3).select();
xlim([0 max(appsizes_mm2)]); ylim([-0.5 ceil(max(appcurr_nA))])
line([0 max(appsizes_mm2)], polyval(pfit, [0 max(appsizes_mm2)]),...
    'Color', 'r')
line(appsizes_mm2, appcurr_nA, 'LineStyle', 'none','Marker', '.',...
    'Color', 'k', 'MarkerSize', 10)
tstr1 = sprintf('Photodiode current per unit area: %2.2f nA/mm^2', pfit(1));
tstr2 = sprintf('Mean irradiance = %G W/mm^2', irradiance_Wpermm2/2);
tstr3 = sprintf('Mean photon flux = %G photons/(s*mm^2)', photon_flux/2);
title({tstr1, tstr2, tstr3});
ylabel('Photodiode current (nA)')
xlabel('Area (mm^2)')

%%
Gid           = fopen(pstr.dpgammatable);
gammatabledat = textscan(Gid,'%f\t%u\t%f\t%s');
fclose(Gid);

gammarange   = mean(reshape(gammatabledat{1}, [256, 2]), 2);
gammacurr_nA = scalePhotodiodeToNanoAmps(gammatabledat{3}, gammatabledat{4});
gammacurr_nA = mean(reshape(gammacurr_nA, [256, 2]), 2);
gammacurve   = gammacurr_nA/max(gammacurr_nA);


% now plot the gammatable
p(1, 4).select();
axis equal; xlim([0 1]); ylim([0 1]);
line(gammarange, gammacurve, 'LineWidth', 2, 'Color', 'k')
line(gammacurve, gammarange, 'LineWidth', 2, 'Color', 'g')
xlabel('Input intensity')
ylabel('Photodiode current (norm.)')
title('Gamma table linearity')
legend({'gamma', 'inverted'}, 'Location', 'southeast')
%%


for ispecies = 1:Nspecies
    
    % calculate isomerizations
    [opsin_absorbance, mouse_opsin_colarea_um2, opsinLabels, opsinCols] = ...
        speciesOpsins(specieslist{ispecies}, wavelengths);

    % integrate photos with respect to the opsin absorbances
    opsin_effective_photons_per_smm2 = opsin_absorbance * photons_per_smm2_per_nm';

    photoiso_rate = opsin_effective_photons_per_smm2 .* mouse_opsin_colarea_um2'* 1e-6;
    
    p(2, ispecies).select();
    pbaspect([2 1 1])
    xlim([start_nm end_nm]);
    ylim([0 1])
    for iopsin = 1:size(opsin_absorbance,1)
        line(wavelengths, opsin_absorbance(iopsin, :), ...
            'Color', opsinCols{iopsin}, 'LineWidth', 2)
    end
    xlabel('Wavelength (nm)')
    ylabel('Relative absorbance')
    legend(opsinLabels,'Location', 'northeast')
    
    opsinstr = cell(numel(opsinLabels), 1);
    for iopsin = 1:size(opsin_absorbance,1)
        opsinstr{iopsin} = sprintf('%s: %d, ', ...
            opsinLabels{iopsin}, round(photoiso_rate(iopsin)/2));
    end
    
    txtrates = strjoin(opsinstr);
    txtrates(end-1:end) = [];
    tstr = sprintf('%s\nmean iso rates\n%s', ...
        specieslist{ispecies}, txtrates);
    
    title(tstr);
end

rez.irradiance_Wpermm2 = irradiance_Wpermm2/2;
rez.photon_flux        = photon_flux/2;

p.export(fullfile(pstr.mainpath, 'res.png'), sprintf('-w%d',fw*10),sprintf('-h%d',fh*10), '-rp');

end

function pathstruct = getPathsFromUser()

[files, path, indx] = uigetfile( ...
{'*.txt;*.csv', 'all relevant files (*.txt,*.csv)';
   '*.csv','csv files (*.csv)'; ...
   '*.txt','text files (*.txt)';
   '*.*',  'all Files (*.*)'}, ...
   'Select calibration files', 'Multiselect', 'on');
%--------------------------------------------------------------------------
assert(numel(files) == 4);
pathstruct.mainpath = path;
%--------------------------------------------------------------------------
% get spectrometer paths
csvinds = find(contains(files, 'csv'));
% get blank apperture
blankind = contains(files(csvinds), {'blank', 'black'});
% set paths
pathstruct.dpspectrumblank = fullfile(path, files{csvinds(blankind)});
pathstruct.dpspectrum      = fullfile(path, files{csvinds(~blankind)});
%--------------------------------------------------------------------------
% get remaining paths
txtinds = find(contains(files, 'txt'));
% get gammatable
gammaind = contains(files(txtinds), {'gamma', 'gama'});
% set paths
pathstruct.dpgammatable    = fullfile(path, files{txtinds(gammaind)});
pathstruct.dpaperturetable = fullfile(path, files{txtinds(~gammaind)});
%--------------------------------------------------------------------------
% get diode file inside folder
pathstruct.dpdiodesens = 'FitSensitivityPhotodiode200to1100nm.txt';
%pathstruct.dpdiodesens = 'pda100a_data.mat';
%--------------------------------------------------------------------------
end



