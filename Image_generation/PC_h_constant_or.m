%Author: Hermine Biermé

% Definition of an AFBF, for which
% h is constant on the conic section (or-alp,or+alp) with 0<alp<=pi/2
% and with a orientation parameter -pi/2<or<=pi/2
% and null outside.

function [ang,c,h] = PC_h_constant_or(H,alp,or)
    if alp <= 0 || alp > pi/2
        fprintf('Error: alp out of range. 0<alp<=pi/2.'); return;
    end
    h=H; c=1;
    ang=[or-alp,or+alp];
    if alp<pi/2
        if or + alp < pi/2 && or - alp > -pi/2
            ang=[-pi/2, ang, pi/2]; c=[0 c 0]; h=[H h H]; 
        elseif or + alp == pi/2
            ang=[-pi/2, ang]; c=[0 c]; h=[H h]; 
        elseif or - alp == -pi/2
            ang=[ang, pi/2]; c=[c 0]; h=[h H];
        else
            ang=[-pi/2, or+alp-pi, or-alp, pi/2]; c=[c 0 c]; h=[h H h]; 
        end
    else
        ang=[-pi/2, pi/2];
    end
