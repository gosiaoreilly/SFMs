% Malgorzata O'Reilly 2023.
% See the text file "instructions_and_conditions_of_use"
% for the conditions of use and how to use.

% We compute the probability density psi(t)_{ij} that
% the SFM first returns to level 0 at time t and does in phase j,
% given start from level 0 in phase i.

% We apply the LST inversion algorithm by Den Iseger in 2006
% and the code in Toutain et al. in 2011,
% adapted to the analysis of the SFMs here.

% References:

% Den Iseger, P., 2006, “Numerical Transform Inversion Using Gaussian
% Quadrature,” Probability in the Engineering and Informational Sciences, 20, pp. 1–44.

% Numerical inversion of Laplace transform for time resolved thermal characterization experiment
% Toutain J., Battaglia J.-L., Pradere C., Pailhes J., Kusiak A., Aregba W., Batsale J.-C.
% (2011) Journal of Heat Transfer, 133 (4) , art. no. 044504.

function [ft]=DenIseger_scalar(t,phasei,phasej)
i=sqrt(-1);
% Lf: anonymous function for the Laplace transform calculation
% t: column vector of times from 0 to tmax
% ft: inverse transform same size as t
nt=numel(t); dt=max(t)/(nt-1);
% parameter settings
N=8*nt; a=44/N;
% Numerical values of lambda and beta see Ref7 p29
li=[0, 6.28318530717958, 12.5663706962589,18.8502914166954, 25.2872172156717, 34.2969716635260, 56.1725527716607, 170.533131190126];
bi=[1, 1.00000000000004, 1.00000015116847, 1.00081841700481, 1.09580332705189, 2.00687652338724, 5.94277512934943, 54.9537264520382];
c=2*i*pi/N; li=a+i.*li; k=0:N; %i=sqrt(?1)!
% the following line may be replaced by a loop on n=size(li) in order
% to avoid memory failure for huge values of nt=size(t)approx 10^6
[k,li]=meshgrid(k,li);
s=(li+c*k)/ dt;
ft=real( Lf_scalar(s,phasei,phasej) );
ft=4*bi*ft/dt;
ft(1)=0.5*(ft(1)+ft(N+1));
% discrete Fourier’s inversion
ft=ifft((ft(1:N)));
ft=real(exp(a.*(0:nt-1)).*(ft(1:nt)));
return