% acceptorFLIMFRET_batch
% Weiyue Chen, April 2013
% Knowing the phasors for donor channel and FRET channel, calculate the
% phasors correspond to acceptor fluorescence only.
function [a,b,G_AM,S_AM] = acceptorFLIMFRET_batch(G_dd,S_dd,G_da,S_da,Ain,Anon)
%Ain:active acceptor
%Anon:passive acceptor
%AM:The donor excitation acceptor emission phasor cloud position if without
%   bleedthrough (phasors corresponding to fluorescence only from acceptor.)
%y=ax+b is the line connecting DD and DM
%For AM,DM,DA,Ain,Anon: [g coordinate; s coordinate]
interAA = det([Ain(1),Ain(2);Anon(1),Anon(2)]);
interDMDA = G_dd.*S_da-S_dd.*G_da;
denomin = (Ain(1)-Anon(1)).*(S_dd-S_da)-(Ain(2)-Anon(2)).*(G_dd-G_da);
G_AM = (interAA.*(G_dd-G_da)-interDMDA.*(Ain(1)-Anon(1)))./denomin;
S_AM = (interAA.*(S_dd-S_da)-interDMDA.*(Ain(2)-Anon(2)))./denomin;
a = (S_dd-S_da)./(G_dd-G_da);
b = interDMDA./(G_dd-G_da);
end