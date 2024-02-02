% aFretSingle
% WeiYue Chen, Aug 2013
% Knowing the active donor lifetime and the passive acceptor lifetime,
% calculate the phasor position of the active acceptor.
function [gAfret,sAfret]=aFretSingle(Dfret,tauA,w)
% Dfret is the position of interacting donor
% A1 and A2 is the position of the non-interacting acceptor components
% w is the laser repetition frequency.
gA=1/((w*tauA)^2+1);
sA=(w*tauA)/((w*tauA)^2+1);
tauD=Dfret(2)./(Dfret(1).*w);
gAfret=(gA*tauA-Dfret(1)*tauD)./(tauA-tauD);
sAfret=(sA*tauA-Dfret(2)*tauD)./(tauA-tauD);
end