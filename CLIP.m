%************************************CLIP.m********************************
%   reduces an angle to within -pi and +pi
%   copyright (c) D. Kajfez 2010
function y=CLIP(x)
while x > pi
    x=x-2*pi;
end
while x < -pi
    x=x+2*pi;
end
y=x;
%disp(['clip.m called, angle= ' num2str(y*190/pi) ' deg']);