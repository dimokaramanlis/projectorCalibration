

function T_absorbance = opsinSpecta(wavelen, opsinpeaks, nomogramsource, varargin)
%
%%% opsinSpecta %%%
%
%
% This function compute normalized absorbance according to various nomograms
% for different opsins at wide range of wavelength. The nomogram options are
% the standard ones use in the retina field.
%
%================================Inputs====================================
%
%   wavelen : input wavelength in nm (e.g. from 100 1000).
%   opsinpeaks : peak absorptance sensitivity of an opsin molecule (e.g.for
%                mouse opsins are [360 496 508] for S-Cone, M-Cone and Rods.
%   nomogramsource : list of nomograms :
%                                          Baylor
%                                          Dawis
%                                          Govardovskii (Default)
%                                          Lamb
%                                          StockmanSharpe
%
%================================Output====================================
%
%   T_absorbance : ormalized absorbance according to selected nomogram.
%
% written by Mohammad, 05.07.2019 based on the PhotopigmentNomogram
% function from Psychtoolbox-3. The sub-functions are slightly modified
% version of the Psychtoolbox-3 functions written by dhb 7/11/03. For more
% info check http://psychtoolbox.org

if (nargin < 3 || isempty(nomogramsource)),	nomogramsource = 'Govardovskii';            end
if length(opsinpeaks)> 1 && isrow(opsinpeaks), opsinpeaks = transpose(opsinpeaks);      end

switch lower(nomogramsource)
    
    case {'govardovskii','g','govar'}
        vislim = false; % set to true to limit the output to visual range
        T_absorbance = Govardovskii_Nomogram_local(wavelen, opsinpeaks, vislim);
        
    case {'dawis','d','daw'}
        T_absorbance = Dawis_Nomogram_local(wavelen, opsinpeaks);
        
    case {'baylor','b','bay'}
        T_absorbance = Baylor_Nomogram_local(wavelen, opsinpeaks);
        
    case {'lamb','l', 'lam'}
        T_absorbance = Lamb_Nomogram_local(wavelen, opsinpeaks);
        
    case {'stockmansharpe','stockman_sharpe','stockman sharpe','stoc','s','ss','stocsharp','sharp'}
        T_absorbance = StockmanSharpe_Nomogram_local(wavelen, opsinpeaks);
        
    otherwise
        error('Da Fuck is this shit! :Unknown source for photopigment nomogram');
end


end

%--------------------------------------------------------------------------------------------------%
%----------                               sub-functions                                  ----------%
%--------------------------------------------------------------------------------------------------%

function T_absr = Govardovskii_Nomogram_local(wls,lambdaMax, varargin)
% T_absorbance = GovardovskiiNomogram(S,lambdaMax)
%
% Compute normalized absorbance according to the
% nomogram provided in Victor I. Govardovskii et al., 2000,
% Visual Neuroscience, Vol. 17, pp. 509-528.
%
% The polynomial approximation has two bands.
% Alpha-band, the major band of A1 pigments is taken from
% Equation (1) and equation (2) in wavelength range [350,700] nm
% (see page 515 in paper)
%
% Beta-band, the minor band of A1 pigments is taken from
% Equation (3) and equation (5a, 5b) in wavelgnth range [350,700] nm
% (see page 516 in paper)
%
% T_absorbance contains the absorbance.
% This is specified according to the CVRL web page
% definition absorbance = log(I_incident/I_transmitted),
% so the numbers are all positive.  The peak is normalized
% to one.
%
% For low density, the absorbance has the same shape as
% the spectral senstivity.
%
% The result is in quantal units, in the sense that to compute
% absorptions you want to incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta().
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% 03/08/2002 ly  Wrote it starting from DawisNomogram.

% Valid range of wavelengh for A1-based visual pigments
% (see page 516 in paper)

% limit the output to visual range
if nargin > 2, vislim = logical(varargin{1});    else, vislim = false; end

Lmin = 330;
Lmax = 700;

% Valid range of lambdamax value
% (taken from figure 4, page 514 in paper)
% this is commented out for convience but used in the paper and psych
% toolbox, it has no effect on the output.
% to activate this, comment out the lines below and then comment out the if
% else condition in the main for loop.

% lmaxLow = 350;
% lmaxHigh = 600;

% Alpha-band parameters
% A = 69.7, B = 28, C = 14.9, D = 0.674 for equation (1) (page 515)
A_B_C = [69.7, 28, -14.9];
D = 0.674;

% b = 0.922, c = 1.104 for equation (1) (page 515)
b_c = [0.922, 1.104];

% Beta-band parameters
% Abeta = 0.26 for equation (4) (page 516)
Abeta = 0.26;

% Get wls argument.
%wls = calib.MakeItWls(S);
if ~iscolumn(wls), wls = transpose(wls); end    % for multiplication it has to column

[nWls,nil] = size(wls);
[nT,nil] = size(lambdaMax);
T_absr = zeros(nT,nWls);

for ii = 1:nT
    theMax = lambdaMax(ii);
    % if (theMax > lmaxLow && theMax < lmaxHigh)
    
    % alpha-band polynomial
    %
    % Parameter a depends on lambdamax, see equation (2) (page 515)
    a = 0.8795 + 0.0459*exp(-(theMax-300)^2/11940);
    a_b_c = [a, b_c];
    
    x = theMax./wls;
    
    % midStepN, N = 1, 2, ... are the middle steps in the caculation.
    midStep1 = exp (ones(nWls,1)*(A_B_C.*a_b_c) - x*A_B_C);
    midStep2 = sum(midStep1,2) + D;
    
    % Result of equation (1) (page 515)
    S_x = 1./midStep2;
    
    % Beta-band polynomial
    
    % Parameter bbeta depends on lambdamax, see equation (5b) (page 516)
    bbeta = -40.5 + 0.195*theMax;
    
    % Conversion of lambdamax to parameter lambdaMaxbeta, see equation (5a) (page 516)
    lambdaMaxbeta = 189 + 0.315*theMax;
    
    % midStepN, N = 1, 2, ... are the middle steps in the caculation.
    midStep1 = -((wls - lambdaMaxbeta * ones (nWls,1)) / bbeta).^2;
    midStep2 = Abeta * exp (midStep1);
    
    % Result of equation (4) (page 516)
    S_beta = midStep2;
    
    % alpha band and beta band together.
    T_absr(ii,:) = (S_x + S_beta)';
    
    if vislim
        % Zero sensitivity outsize valid range.
        index = find(wls < Lmin);
        T_absr(ii,index) = zeros(size(index))';
        index = find(wls > Lmax);
        T_absr(ii,index) = zeros(size(index))';
    end
    % else
    %    error('Lambda Max %g not in range of nomogram\n',theMax);
    % end
    
end
end

%--------------------------------------------------------------------------------------------------%

function T_absr = Dawis_Nomogram_local(wls,lambdaMax)
% [T_absorbance] = DawisNomogram(S,lambdaMax)
%
% Compute normalized absorbance according to the
% nomogram provided in Dawis, 1981, Vision Research,
% Vol. 21, pp. 1427-1430.
%
% T_absorbance contains the absorbance.
% This is specified according to the CVRL web page
% definition absorbance = log(I_incident/I_transmitted),
% so the numbers are all positive.  The peak is normalized
% to one.
%
% For low density, the absorbance has the same shape as
% the spectral senstivity.
%
% The nomogram is shifted along the wavelength axis
% using a multiplicative rather than additive procedure.
%
% The result is in quantal units, in the sense that to compute
% absorptions you want to incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta().
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% 10/30/97 dhb  Wrote it.
% 07/01/03 dhb  Add computation of T_absorbance.

% These are the coefficients for the polynomial
% approximation, taken from Table 1.  We implement
% the A1-based pigment nomogram, not the A2.
Lmax = [432 ; 502 ; 562];
bN = [0.346325  -35.1001  63.3807 -125.466  10962.7  -16244.8 -210671.0   -23776.9 ; ...
    -0.0106836 -28.28   148.133  -498.627  -1457.94  12799.4    -789.371 -60749.2 ; ...
    -0.228262  -22.9974  87.0027 -636.336   2624.45   4948.6  -45944.6    65688.5 ];
lmaxLow =  [410 ; 470 ; 530];
lmaxHigh = [470 ; 530 ; 610];
L1 = [380 ; 400 ; 430 ];
L2 = [510 ; 620 ; 690 ];

% Get wls argument.
%wls = calib.MakeItWls(S);

[nWls,nil] = size(wls);
[nT,nil] = size(lambdaMax);
T_absr = zeros(nT,nWls);
wlsum = wls/1000;

for i = 1:nT
    theMax = lambdaMax(i);
    if (theMax >= lmaxLow(1) && theMax <= lmaxHigh(1))
        whichval = 1;
    elseif (theMax > lmaxLow(2) && theMax <= lmaxHigh(2))
        whichval = 2;
    elseif (theMax > lmaxLow(3) && theMax <= lmaxHigh(3))
        whichval = 3;
    else
        error('Lambda Max %g not in range of nomogram, acceptable range is %g to %g \n',...
            theMax, lmaxLow(whichval), lmaxHigh(whichval));
    end
    wlsVec = (theMax ./ wls') - 1;
    logS = zeros(1,nWls);
    for k = 1:8
        logS = logS + bN(whichval,k)*wlsVec.^k;
    end
    T_absr(i,:) = logS;
    
    % Zero sensitivity outsize valid range.  I shift the
    % range in the table to slide multiplicatively with
    % theMax.
    zeroLow = theMax*L1(whichval)/Lmax(whichval);
    zeroHigh = theMax*L2(whichval)/Lmax(whichval);
    T_absr(i,wls < zeroLow) = -Inf;
    T_absr(i,wls > zeroHigh) = -Inf;
end
T_absr = 10.^T_absr;

end

%--------------------------------------------------------------------------------------------------%

function T_absr = Baylor_Nomogram_local(wls,lambdaMax)
% T = BaylorNomogram(S,lambdaMax)
%
% Compute spectral sensitivities according to the
% nomogram provided in Baylor, Nunn, and Schnapf, 1987.
%
% The result is in quantal units, in the sense that to compute
% absorptions you want to incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta().
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% 6/22/96  dhb  Wrote it.
% 10/16/97 dhb  Add comment about energy units.

% These are the coefficients for the polynomial
% approximation.
aN = [-5.2734 -87.403 1228.4 -3346.3 -5070.3 30881 -31607];

% Get wls argument.
%wls = calib.MakeItWls(S);

[nWls,nil] = size(wls);
[nT,nil] = size(lambdaMax);
T_absr = zeros(nT,nWls);
wlsum = wls/1000;

if max(wls) > 800
    error('Maximum wavelength range %g not in range of nomogram, acceptable range is %g to %g \n',...
        max(wls), 100,800);
end

for i = 1:nT
    if lambdaMax(i) < 420
        error('Lambda Max %g not in range of nomogram, acceptable range is %g to %g \n',...
            lambdaMax(i), 420,900);
    end
    wlsVec = log10( (1 ./ wlsum)*lambdaMax(i)/561)';
    logS = aN(1) + aN(2)*wlsVec + aN(3)*wlsVec.^2 + aN(4)*wlsVec.^3 + ...
        aN(5)*wlsVec.^4 + aN(6)*wlsVec.^5 + aN(7)*wlsVec.^6;
    T_absr(i,:) = 10.^logS;
end

end

%--------------------------------------------------------------------------------------------------%

function T_absr = Lamb_Nomogram_local(wls,lambdaMax)
% T = LambNomogram(S,lambdaMax)
%
% Compute spectral sensitivities according to the
% nomogram provided in Lamb, 1995, Vision Research,
% Vol. 35, pp. 3083-3091, equation 2'.
%
% The result is in quantal units, in the sense that to compute
% absorptions you want to incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta().
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% 8/20/98 dhb  Wrote it.

% These are the coefficients for Equation 2.
a = 70; b = 28.5; c = -14.1;
A = 0.880; B = 0.924; C = 1.104; D = 0.655;

% Get wls argument.
%wls = calib.MakeItWls(S);
[nWls,nil] = size(wls);
[nT,nil] = size(lambdaMax);
T_absr = zeros(nT,nWls);

for i = 1:nT
    theMax = lambdaMax(i);
    wlarg = theMax ./ wls';
    T_absr(i,:) = 1 ./ ( exp( a*(A-wlarg) ) + exp( b*(B-wlarg) ) + exp( c*(C-wlarg) ) + D );
    T_absr(i,:) = T_absr(i,:)/max(T_absr(i,:));
end

end

%--------------------------------------------------------------------------------------------------%

function T_absr = StockmanSharpe_Nomogram_local(wls,lambdaMax)
% T_absorbance = StockmanSharpeNomogram(S,lambdaMax)
%
% Compute normalized absorbance according to the
% nomogram provided by Stockman and Sharpe:
%   See Stockman & Sharpe (2000), p. 1730 or http://www.cvrl.org/.
%
% The lmax values of the fitted templates that best fits
% to the Stockman and Sharpe (2000) S-, M- and L-cone photopigment
% spectra are 420.7, 530.3 and 558.9 nm for the
% S-, M- and L-cones, respectively;
%
% The result is in quantal units, in the sense that to compute
% absorptions you want to incident spectra in quanta.
% To get sensitivity in energy units, apply EnergyToQuanta()
% (not QuantaToEnergy, because here you are converting sensitivity
% rather than spectra.)
%
% Argument lambdaMax may be a column vector of wavelengths.
%
% This routine converts the log10 absorbance computed by the nomogram
% formulae into absorbance.
%
% Note from DHB. By eye, this nomogram gives a good fit to the log10 LMS pigment
% absorbance spectra when you plot them on a log scale over 8 log units,
% as in Stockman & Sharpe (2000), Figure 12.  The fit does not look quite so
% good in my hands on a linear scale or an expanded log scale.
% I think the deviations occur in part because the tabulated
% photopigment absorbances for the L are a mixture of the ser/ala
% pigments, and the nomogram was developed in part on the basis of fitting these
% as if they corresponded to a single lambda max value. It may also be that
% the template was built by minimizing the error in log sensitivity, and this
% would more heavily weight the long wavelength limbs of the pigments, where
% the linear sensitivity is essentially zero.  In any case, though, if you
% are working in the land of the CIE 170-1:2006 fundamentals, this is the
% probably the best current nomogram to use.
%
% See ComputeCIEConeFundamentals, CIEConeFundamentalsTest,
% FitConeFundamentalsFromNomogram, FitConeFundamentalsTest
%
% 5/8/99	dhb  Started writing it.
% 10/27/99	dhb  Added error return to prevent premature use of this routine.
% 7/18/03   dhb  Finished it off.
% 8/13/11   dhb  Improved comments.  Double check polynomial coefficients.

% Parameters
a = -188862.970810906644;
b = 90228.966712600282;
c = -2483.531554344362;
d = -6675.007923501414;
e = 1813.525992411163;
f = -215.177888526334;
g = 12.487558618387;
h = -0.289541500599;

% Set up and apply formula
%wls = calib.MakeItWls(S)';
nWls = length(wls);
nT = length(lambdaMax);
T_absr = zeros(nT,nWls);
for i = 1:nT
    % Get lambda max
    theMax = lambdaMax(i);
    
    % Need to normalize wavelengths
    logWlsNorm = log10(wls)-log10(theMax/558);
    
    % Compute log optical density
    logDensity = a + b*logWlsNorm.^2 + c*logWlsNorm.^4 + d*logWlsNorm.^6 + ...
        e*logWlsNorm.^8 + f*logWlsNorm.^10 + g*logWlsNorm.^12 + h*logWlsNorm.^14;
    %logDensity = logDensity;
    T_absr(i,:) = 10.^logDensity;
end

end
