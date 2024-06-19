%Author: Hermine Biermé

% Definition of an AFBF, using the turning bands method and Frederic Richard's
% matlab code (see PC_TurningBandsV4.m)
% for which h is constant on the conic section (-alp,alp) with 0<alp<=pi/2
% and null outside.
function [ang,c,h] = PC_h_constant(H,alp)
h=H; c=1;
ang=[-alp,alp];
if alp<pi/2, ang=[-pi/2, ang, pi/2]; c=[0 c 0]; h=[H h H]; end
