% AUTHOR: R. Soulard, XLIM Lab. CNRS 7252 University of Poitiers, France
%         raphael.soulard@univ-poitiers.fr
%
% 2015 
%
function b = anglewrap(a,opt)
% wraps an angle data 'a' into some restricted interval 
% by adaptively adding or substracting %pi or 2*%pi
% usage : b=anglewrap(a,'[-pi/2;pi/2]');
pipi = a - 2*pi*round(a/(2*pi)); % set a in [-pi;pi]
switch opt,
  case '[-pi;pi]', b = pipi;
  case '[0;2pi]', b=pipi; b(b<0)=b(b<0)+2*pi;
  case '[-pi/2;pi/2]', b=pipi;
    b(b<-pi/2)=b(b<-pi/2)+pi;
    b(b>+pi/2)=b(b>+pi/2)-pi;
  case '[0;pi]', b=pipi; b(b<0)=b(b<0)+pi;
  otherwise b=a; fprintf('anglewrap(): bad argument!');
end
