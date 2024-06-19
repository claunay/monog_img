%
% This source code is freely distributed from the "colormonogenic" website:
% http://xlim-sic.labo.univ-poitiers.fr/projets/colormonogenic/
% published in 2018,
% which presents the main research results by 
% Raphaël Soulard & Philippe Carré,
% from the XLIM Laboratory (UMR CNRS 7252),
% University of Poitiers, France.
%
% Author: R. Soulard.
%

% Subfunction for the radial filterbank 'mwt_radial.m'
% Compute Fourier domain responses of low-pass and high-pass radial filters
% (No dependencies)
%

% Supported wavelet types 'typ':
% {'GaussianHP'}
% {'GaussianLP'}
% {'Storath',order} (order in {0..5})
% {'UnserSimonc'}
% {'UnserMeyer'}
% {'UnserPapadakis'}
% {'UnserShannon'}
% {'PadU1D'}

function [LP,HP] = mwt_get_filters( RHO , typ )
%
% RHO is a 2d matrix containing the Euclidean norm of Fourier coordinates,
% which takes values in [0;pi*sqrt(2)].
% The computed Fourier response is radial (isotropic) and in [0;1]
%
switch typ{1},
  
  case 'GaussianHP', % -------------------------------------------------- %
    HP  =  1 - exp( -(RHO.^2)/2 ); % Gaussian-based high-pass filter
    LP  =  sqrt( 1 - HP.^2 );               % Complementary low-pass
  case 'GaussianLP', % -------------------------------------------------- %
    LP  =  exp( -(RHO.^2)/8 );       % Gaussian low-pass filter
    HP  =  sqrt( 1 - LP.^2 );          % Complementary high-pass filter
  case 'Storath', % ----------------------------------------------------- %
    % Monogenic filterbank from Held et al.
    % The sign of the Riesz transform will be different from original
    % definition (up to a minus sign, for consistency).
    % The original reconstruction method (from primary part only) will also
    % differ from the reference.
    if length(typ)>1,
      w_order = typ{2};
    else
      w_order = 2; % parameter 'order' in {0,1,2,3,4,5}
    end
    
    % Evaluate the polynomial 'q':
    x = 4*RHO/pi;
    switch w_order
    case 0, q = 1/2 - 1/4 *x;
    case 1, q = 1/2 - 1/4 *x;
    case 2, q = 8 - 30 *x + 45 *x.^2 - 65/2 *x.^3 + 45/4 *x.^4 - 3/2 *x.^5;
    case 3, q = -52 + 280 *x - 630 *x.^2 + 770 *x.^3 - 2205/4 *x.^4 ...
                + 231 *x.^5 - 105/2 *x.^6 + 5 *x.^7;
    case 4, q = 368 - 2520 *x + 7560 *x.^2 - 13020 *x.^3 + 14175 *x.^4 ...
                - 20223/2 *x.^5 + 4725 *x.^6 - 1395 *x.^7 + 945/4 *x.^8 ...
                - 35/2 *x.^9;
    case 5, q = - 2656 + 22176 *x - 83160 *x.^2 + 184800 *x.^3 ...
                - 270270 *x.^4 + 273042 *x.^5 - 388773/2 *x.^6 ...
                + 97515 *x.^7 - 135135/4 *x.^8 + 7700 *x.^9 ...
                - 2079/2 *x.^(10) + 63 *x.^(11);
    otherwise
      error('Choose order less than 6')
    end
    % NB: When RHO rises from pi/4 to pi/2   ,   q falls from 1/4 to 0
    %     (and then cos( 2*pi*q ) rises from 0 to 1 )
    % Higher w_order increases response's steepness
    HP = cos( 2*pi*q ); % High-pass: rising transition between pi/4 and pi/2
    HP( RHO < pi/4 ) = 0;
    HP( RHO > pi/2 ) = 1;
    LP = sqrt( 1 - HP.^2 ); % Low-pass: falling transition between pi/4 and pi/2
  case 'UnserSimonc', % ------------------------------------------------- %
    HP = cos( (pi/2)*log2(2*RHO/pi) ); % High-pass rising from 0 to 1 in [pi/4;pi/2]
    HP( RHO < pi/4 ) = 0;
    HP( RHO > pi/2 ) = 1;
    LP = sqrt( 1 - HP.^2 );
  case 'UnserMeyer', % -------------------------------------------------- %
    t = 4*RHO/pi-1;
    q = (t>0).*(t<=1).*(t.^4 .*(35 - 84.*t + 70*t.^2 -20*t.^3));
    LP = cos(  (pi/2)*q   );
    LP( RHO < pi/4 ) = 1;
    LP( RHO > pi/2 ) = 0;
    HP = sqrt( 1 - LP.^2 );
  case 'UnserPapadakis', % ---------------------------------------------- %
    LP = sqrt(  0.5+0.5*cos(  5*RHO - 3*pi/2  ) ) ;
    LP( RHO < 3*pi/10 ) = 1;
    LP( RHO > pi/2 ) = 0;
    HP = sqrt(1- LP.^2);
  case 'UnserShannon', % ------------------------------------------------ %
    LP = (RHO <  (pi/2));
    HP = (RHO >= (pi/2));
  case 'PadU1D', % ------------------------------------------------------ %
    LP = RHO*0; % Init lowpass
    idx = RHO>=pi/4 & RHO<=pi/2;
    LP(idx) = sqrt((  (log2(pi./(2*RHO(idx))) + 0.6).^4 - 0.1296  )/6.424);
    LP(RHO<pi/4) = 1;
    HP = sqrt( 1 - LP.^2 );
end

