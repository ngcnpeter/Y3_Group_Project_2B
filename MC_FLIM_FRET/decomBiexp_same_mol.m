% decomBiexp_same_mol
% WeiYue Chen, Aug 2013
% For a bi-exponential decay phasor, assume that the MOLECULE NUMBER
% of the two lifetime components are the same, and do the decomposition of
% the phasor and obtain the two single exponential decay components.
function [A,B,rt1,rt2,tau1,tau2]=decomBiexp(g,s,w)
% g is the x-axis coordinate of the mixture phasor.
% s is the y-axis coordinate of the mixture phasor.
% w is the laser frequency (Fourier transform frequency)
% The line connecting the two roots: y=ax+b
% rt1 is the coordinate of the first root, rt1(1)is the x coordinate, 
% rt1(2) is the y coordinate. The same for rt2.
% tau1 and tau2 is the corresponding lifetime of the rt1 and rt2
syms a b root1 root2 Tau1 Tau2
b=s-a*g;
root1=(1-2*a*b+sqrt(1-4*a*b-4*b^2))/(2*(a^2+1));
Tau1=(a*root1+b)./(root1*w);
root2=(1-2*a*b-sqrt(1-4*a*b-4*b^2))/(2*(a^2+1));
Tau2=(a*root2+b)./(root2*w);
A=solve([char(Tau1),'*(',char(root1),'-',num2str(g),')=(',num2str(g),'-',char(root2),')*',char(Tau2)],'a');
A=double(A);
B=s-A*g;

rt1(1)=(1-2*A*B+sqrt(1-4*A*B-4*B^2))/(2*(A^2+1));
rt1(2)=A*rt1(1)+B;
rt2(1)=(1-2*A*B-sqrt(1-4*A*B-4*B^2))/(2*(A^2+1));
rt2(2)=A*rt2(1)+B;

tau1=rt1(2)./(rt1(1).*w);
tau2=rt2(2)./(rt2(1).*w);
end