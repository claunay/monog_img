% AUTHOR: R. Soulard, XLIM Lab. CNRS 7252 University of Poitiers, France
%         raphael.soulard@univ-poitiers.fr
%
% 2015 
%
function [RHO,RZ]=FFT_radial(siz)
% Creates an MxN grid of Fourier coordinates 
% - First sample is the "DC coordinate" (no 'fftshift' done)
% - Frequencies span [-pi;pi]^2

% 'siz' : size in the form '[M,N]'
% 'RHO' : absolute frequency or radial frequency (distance from DC)
% 'RZ'  : 'orientation', embedded in a complex exponential
% 'RZ' turns out to coincide with the Riesz-transform frequency response.
% The singular value of riez at DC is set to 1 for convenience
% (should be 0 in theory).
M = siz(1);
N = siz(2);
[W1,W2] = ndgrid(  (-floor(M/2):ceil(M/2)-1)/M  ,...
                   (-floor(N/2):ceil(N/2)-1)/N  );
RHO  =  ifftshift(        2*pi*sqrt( W1.^2 + W2.^2 )  ); % Radial freq. coord.
RZ = ifftshift(  -1i * exp(  1i*atan2(W2,W1)  )    ); % Riesz filter
