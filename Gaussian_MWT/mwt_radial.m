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

% The function performs a greyscale Monogenic Wavelet Transform
% based on radial filters defined in the FFT domain.
% Both undecimated and pyramid settings are possible, 
% but note that the pyramid setting will not necessarily provide perfect
% reconstruction (depending on the type of wavelet).
%
% N.B: The pyramid setting is still in progress... 
%      (we have been only using the undecimated setting for our work)
%      (Note that Gaussian-based wavelets do not fulfill the tight-frame
%       property, so are not eligible for reversible pyramid filterbank)
%
% Several types of wavelets are covered:
% - Monogenic wavelets from [Held, Storath et al. "Steerable wavelet frames..." IEEE TIP 2010]
% - Steerable wavelets from [Unser & Chenouard "A unifying framework..." SIAM J.Im.Sci. 2013]
% - Alternative Gaussian-based wavelets from 
%  [Soulard & Carre, "Characterization of color images..." IEEE TPAMI 2017]

% Supported wavelet types 'typ':
% {'GaussianHP'}                  [Soulard & Carre IEEE TPAMI 2017]
% {'GaussianLP'}                   Soulard & Carré (unpublished)
% {'Storath',ord} (ord in {0..5}) [Held, Storath et al. IEEE TIP 2010]
% {'UnserSimonc'}                 [Unser & Chenouard SIAM J.Im.Sci. 2013]
% {'UnserMeyer'}                  [Unser & Chenouard SIAM J.Im.Sci. 2013]
% {'UnserPapadakis'}              [Unser & Chenouard SIAM J.Im.Sci. 2013]
% {'UnserShannon'}                [Unser & Chenouard SIAM J.Im.Sci. 2013]
% {'PadU1D'}                      [Unser & Chenouard SIAM J.Im.Sci. 2013]

%
% Dependencies:
%
%   'FFT_radial.m' (no dependency)  
%     % Get Fourier coordinates:
%     [RHO,RZ] = FFT_radial(siz);
%
%   'mwt_get_filters.m' (no dependency)
%     % Get Fourier responses of different wavelets:
%     [LP,HP] = mwt_get_filters( RHO , typ );
%
%   'mwt_get_norm_csts.m'  (no dependency)
%     % Get normalization constants associated to chosen wavelets:
%     norm_csts = mwt_get_norm_csts( typ );
%


function out = mwt_radial( opt , in , L , param )
% Get input parameters:
typ      = param.typ;
norma    = param.norma;
sampling = param.sampling;
noRzHF   = param.noRzHF;
% opt : 'a' for analysis, 's' for 'synthesis'
% in  : MxN 2d matrix input image or (L+1,3)-cell-array input MWT
% L   : integer - number of scales
% typ : string withing a cell (such as {'GaussianHP'}) - type of wavelet
% norma : 'true': subbands are normalized by the energy of equiv. filters
%         'false' : equiv. filters have freq. resp. in [0;1]
% sampling : 'undec' for the undecimated setting
%            'pyramid' for the pyramid setting
%           In both cases, the output MWT is a cell array of size (L+1,3).
%           The size of subbands in each cell depends on the sampling option.
% noRzHF : 'true' : the Riesz transform is not performed at 1st scale
%                   to enhance the decay of monogenic amplitude.
%          'false': classical MWT


% Handle normalization of subbands:
if norma, % normalise to have RMS=1 for equivalent filter at all scales
  norm_csts = mwt_get_norm_csts( typ );
  norm_csts = [norm_csts ; norm_csts(end)*ones(L-length(norm_csts),1) ]; % (in case L is too large)
  norm_csts = norm_csts(1:L);
else % cancel normalization to have "frequency responses in [0;1]":
  norm_csts = (2.^([0:L-1]')); % compensate the "lfband*2"
  %norm_csts = ones(L,1); % (for debugging or energy measurements)
end

%------------------------------------------------------------------------ %
if opt=='a', % ----------------- MWT Analysis --------------------------- %
%------------------------------------------------------------------------ %

% 'in'  is a 2d matrix (input greyscale image)
% 'out' will be a (L+1,3) cell array (output MWT)
  
  lfband = in;           % Init low-frequency band
  out = cell(L+1,3);     % Init Wavelet transform
  [RHO,RZ] = FFT_radial(  size(lfband)  ); % FFT coordinates & Riesz Tr. 
  for scal=1:L,       % span scales from fine to coarse
    
    [LP,HP] = mwt_get_filters( RHO , typ ); % filters' FFT
    LF = fft2(lfband);
    lfband = real(ifft2( LF .* LP  ));  % coarser approximation
    HF = LF .* HP / (norm_csts(scal)*sqrt(2));    % normalized high-pass
                                          % The sqrt(2) term corresponds to
                                          % the monogenic frame setting.
    if scal==1 && noRzHF,
      riez = HF * 0; % no Riesz part at finest scale
      HF = HF * sqrt(2); % compensating amplification of the primary part
    else
      riez = ifft2( HF .* RZ );           % Riesz transform
    end
    out{scal,1} = real(ifft2(HF)); % Primary part
    out{scal,2} = real(riez);      % x-Riesz part
    out{scal,3} = imag(riez);      % y-Riesz part
    
    if strcmp(sampling,'pyramid'),
      lfband = lfband(1:2:end,1:2:end);
      [RHO,RZ] = FFT_radial(  size(lfband)  ); % new FFT coordinates
    else % undecimated setting
      RHO = RHO*2; % Dilatation of the filters
      lfband = lfband*2; % normalization due to the undecimated setting
    end
    
    
  end 
  out{L+1,1} = lfband;               % store coarsest approximation
  % End of MWT analysis.
  

  %test
%------------------------------------------------------------------------ %
elseif opt=='s', % ------------- MWT Synthesis -------------------------- %
%------------------------------------------------------------------------ %

% 'in'  is a (L+1,3) cell array (input MWT) 
% 'out' will be a 2d matrix (output greyscale image)
%       (iteratively built in variable 'lfband')

  lfband = in{L+1,1};                      % Get coarse approximation
  [RHO,RZ] = FFT_radial(  size(lfband)  ); % FFT coordinates & Riesz Tr. 
  
  [M,N]=size(lfband);
  if ~strcmp(sampling,'pyramid'), % undecimated setting
    RHO = RHO * (2^L);
  end
  
  for scal=L:-1:1, % span scales from coarse to fine
    
    if strcmp(sampling,'pyramid'), % pyramid case
      M = M*2;   N = N*2; % Update size
      [RHO,RZ] = FFT_radial( [M,N] ); % Update Fourier coordinates and Riesz filter
      % Upsampling of 'lfband':
      tmp = zeros(M,N);
      tmp(1:2:end,1:2:end) = 4*lfband; % 4x 1/4th of the coeffs 
      lfband = tmp;   clear tmp;
    else % undecimated case 
      RHO = RHO/2; %Dilatation of the filters
      lfband = lfband/2;% Normalization due to the undecimated setting
    end
    
    [LP,HP] = mwt_get_filters( RHO , typ ); % filters' FFT
    HP = HP * norm_csts(scal);   % de-normalize the high-pass filter
    LF = fft2(lfband)  .* LP ;   % low-pass filter the approximation
    PRIM = fft2( in{scal,1}           )  .* HP  ;       % High pass
    
    % New approximation 'lfband':
    if scal==1 && noRzHF, % No Riesz transform at first scale:
      lfband = real(ifft2( LF + PRIM ));                
    else %                % General case:
      % Retrive Riesz part:
      RIEZ = fft2( in{scal,2} + 1i*in{scal,3} )  .* HP .* conj(RZ) ; 
      lfband = real(ifft2(  LF + ( PRIM + RIEZ )/sqrt(2)  )); 
      % (the sqrt(2) term corresponds to the monogenic frame setting)
    end
    
  end
  out = lfband; % Reconstructed greyscale image
  % End of MWT reconstruction.

%------------------------------------------------------------------------ %
else % bad 'opt' arg
  fprintf('mwt - bad arg.: ''opt'' must be ''a'' or ''s''.\n'); return;
end



