% aFret
% WeiYue Chen, Aug 2013
% Knowing the active donor lifetime, the two lifetime components of
% the passive acceptor lifetime, calculate the phasor position of
% the active acceptor.
function [gAfret,sAfret]=aFret(Dfret,A1,A2,w)
% Dfret is the position of interacting donor
% A1 and A2 is the position of the non-interacting acceptor components
% w is the laser repetition frequency.
tau1=A1(2)./(A1(1).*w);
tau2=A2(2)./(A2(1).*w);
tauD=Dfret(2)./(Dfret(1).*w);
fracA1=tau2-tauD;
fracA2=tau1-tauD;
gAfret=(fracA1*A1(1)*tau1+fracA2*A2(1)*tau2-(fracA1+fracA2)*Dfret(1)*tauD)./(fracA1*tau1+fracA2*tau2-(fracA1+fracA2)*tauD);
sAfret=(fracA1*A1(2)*tau1+fracA2*A2(2)*tau2-(fracA1+fracA2)*Dfret(2)*tauD)./(fracA1*tau1+fracA2*tau2-(fracA1+fracA2)*tauD);
end