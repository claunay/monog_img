%Author: Hermine Biermé

% Definition of an AFBF, for which
% h is constant on the conic section (or-alp,or+alp) with 0<alp<=pi/2
% and null outside.
function [ang,c,h] = PC_h_constant_or(H,alp,or)
if abs(or) > 1.175
    fprintf('Error: Choose a smaller rotation angle parameter ''or''.'); return;
end
h=H; c=1;
ang=[or-alp,or+alp];
if or+alp<pi/2, ang=[-pi/2, ang, pi/2]; c=[0 c 0]; h=[H h H]; end
