% tauComp
% Calculate the intersection point of a line and the universal circle
% WeiYue Chen, Aug 2013
function [x,y,tau]= tauComp(a,b,w)
% The fitted line is: y=ax+b
% w: laser frequency (Fourier transform frequency when producing phasor plot)
p=[(a^2+1),(2*a*b-1),b^2];
x=roots(p);
y=a.*x+b;
tau=y./(x.*w);
end