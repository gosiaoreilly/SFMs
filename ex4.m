% Malgorzata O'Reilly 2023.
% See the text file "instructions_and_conditions_of_use"
% for the conditions of use and how to use.

% Example 4: parameters and Psi.

clear all
close all

% State space.
s1=18;
s2=2;
s0=0;

% Rates c_i.
cvec=[ones(1,s1) -ones(1,s2)];
C1=diag(cvec(1:s1));
C2=diag(cvec(s1+1:s1+s2));

% Generator T of the Markov chain.
T11=ones(s1,s1)*10-eye(s1)*10;
T22=ones(s2)*10;
T12=ones(s1,s2)*0.001;
T21=ones(s2,s1)*0.001;
T11=T11-diag(sum(T12,2)+sum(T11,2));
T22=T22-diag(sum(T21,2)+sum(T22,2));
T=[T11,T12;T21,T22];

% Calculate recurrence measure "mu".
pivec=[zeros(1,s1+s2) 1]/[T,ones(s1+s2,1)];
mu = pivec*cvec';

save examplepar.mat T11 T22 T12 T21 s1 s2 s0 C1 C2 mu

% Fluid generator Q (in this example, same as T due to +1/-1 rates c_i).
Q11=inv(C1)*(T11);
Q22=inv(-C2)*(T22);
Q12=inv(C1)*(T12);
Q21=inv(-C2)*(T21);
Q=[Q11,Q12;Q21,Q22];

% Compute Psi matrix.
[Psi, iterationsN]=A4_getPsi(Q11,Q12,Q21,Q22);
Psi
% 'sums of Psrows in Psi'
% sum(Psi,2)


% % Compute psi(t) matrix.
% vect=[1:0.01:2];
% for phasei=1:s1
%     for phasej=1:s2
%         [ft]=DenIseger_scalar(vect,phasei,phasej);
%         plot(vect,ft)
%         hold on
%         xlabel('t','FontSize',14)
%         ylabel('\psi(t)_{ij}','FontSize',14)
%     end
% end