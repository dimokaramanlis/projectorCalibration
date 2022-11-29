function [opsinSpectra, opsinCollAreas, opsinLabels, opsinCols] = ...
    speciesOpsins(currspecies, wavelengths)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

supportedspecies = {'mouse_s', 'mouse_g', 'marmoset'};

assert(sum(contains(supportedspecies, currspecies)) == 1, ...
    'Provided species is not supported')


switch currspecies
    case 'mouse_s'
        [opsinSpectra, opsinCollAreas,opsinLabels, opsinCols] = ...
            mouseOpsins(wavelengths, 's');
    case 'mouse_g'
        [opsinSpectra, opsinCollAreas,opsinLabels, opsinCols] = ...
            mouseOpsins(wavelengths, 'g');
    case 'marmoset'
        [opsinSpectra, opsinCollAreas,opsinLabels, opsinCols] = ...
            marmosetOpsins(wavelengths(:));
end



end

function [opsinspectr, opsin_colarea_um2, opsinnames, opsincols] = ...
    mouseOpsins(wavelengths, temptype)

% numbers from literature for mouse S-opsin, M-opsin and rhodopsin
opsin_peak_nm     = [360, 511, 500];
opsin_colarea_um2 = [0.2, 0.2, 0.5]; 

% collecting areas are the most important, here we use what T. Euler is
% using, see https://github.com/eulerlab/open-visual-stimulator/blob/master/calibration_mouse/stimulator_calibration.ipynb
% cone value of 0.2 comes from Nikonov et al., 2006, rod value was measured
% by G. Field and F. Rieke, 2002 to be 0.52

opsinspectr = opsinSpecta(wavelengths, opsin_peak_nm, temptype);
opsinnames  = {'S', 'M', 'Rh'};
opsincols   = {'m', 'g', 'k'};

end


function [opsinspectr, opsin_colarea_um2, opsinnames, opsincols] = ...
    marmosetOpsins(wavelengths)


% for peaks, see Travis et al 1988 (Vision Research), Tovee et al 1992
% (Vision Research)

opsin_peak_nm     = [423,   543,  556,  563, 499];
opsin_colarea_um2 = [0.37, 0.37, 0.37, 0.37,   1]; 

% collecting areas from macaque, cones: (Schnapf et al 1990), rods:
% (Schneeweis and Schnapf, 1995)

opsinspectr = opsinSpecta(wavelengths, opsin_peak_nm, 'lamb');

opsinnames  = {'S', 'M1', 'M2', 'M3', 'Rh'};
opsincols   = {'b', 'g',   'y',  'r',  'k'};

end


