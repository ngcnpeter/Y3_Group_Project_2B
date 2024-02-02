% decombiexp
% WeiYue Chen, Aug 2013
% For a bi-exponential decay phasor, assume that the INTENSITY CONTRIBUTION
% of the two lifetime components are the same, and do the decomposition of
% the phasor and obtain the two single exponential decay components.
function [Af,Bf,RT1,RT2,tau1,tau2]=decomBiexp(g,s,w)
% g is the x-axis coordinate of the mixture phasor.
% s is the y-axis coordinate of the mixture phasor.
% w is the laser frequency (Fourier transform frequency)
% The line connecting the two roots: y=ax+b
% rt1 is the coordinate of the first root, rt1(1)is the x coordinate, 
% rt1(2) is the y coordinate. The same for rt2.
% tau1 and tau2 is the corresponding lifetime of the rt1 and rt2
syms a b root1 root2
b=s-a*g;
root1=(1-2*a*b+sqrt(1-4*a*b-4*b^2))/(2*(a^2+1));
root2=(1-2*a*b-sqrt(1-4*a*b-4*b^2))/(2*(a^2+1));
A=solve([char(root1),'-',num2str(g),'=',num2str(g),'-',char(root2)],'a');
A=double(A);
    
B=s-A.*g;

rt1(:,1)=(1-2.*A.*B+sqrt(1-4.*A.*B-4.*B.^2))./(2*(A.^2+1));
rt1(:,2)=A.*rt1(:,1)+B;
rt2(:,1)=(1-2.*A.*B-sqrt(1-4.*A.*B-4.*B.^2))./(2.*(A.^2+1));
rt2(:,2)=A.*rt2(:,1)+B;
rt1=rt1';
rt2=rt2';


if size(A,1)==2
   if rt1(2,1)>0 && rt2(2,1)>0
       RT1=rt1(:,1);
       RT2=rt2(:,1);
       Af=A(1);
       Bf=B(1);
   else
       RT1=rt1(:,2);
       RT2=rt2(:,2);
       Af=A(2);
       Bf=B(2); 
   end
end
tau1=RT1(2)./(RT1(1).*w);
tau2=RT2(2)./(RT2(1).*w);

end