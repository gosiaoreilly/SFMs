% Malgorzata O'Reilly 2023.
% See the text file "instructions_and_conditions_of_use"
% for the conditions of use and how to use.

function [psin, iterationsN]=A4_getPsi(Q11,Q12,Q21,Q22)

% For the stopping criterion.
Neps=1e12;
max_iter=100;

Q=[Q11,Q12;Q21,Q22];
s1=size(Q11,1);
s2=size(Q22,1);

% Calculate recurrence measure "mu".
cvec=[ones(1,s1) -ones(1,s2)];
pivec=[zeros(1,s1+s2) 1]/[Q,ones(s1+s2,1)];
mu = pivec*cvec';


% Newton's method (A4) %
A=Q11;
B=Q22;
C=-Q12;
% X = sylvester(A,B,C) solves the Sylvester equation A*X + X*B = C.
psin=sylvester(A,B,C);
oldpsin=psin;
errorN=2;
iterationsN=1;
if errorN>1
    g=1;
    while g<max_iter
        A=Q11+psin*Q21;
        B=Q22+Q21*psin;
        C=-Q12+psin*Q21*psin;
        psin=sylvester(A,B,C);
        iterationsN=iterationsN+1;
        if mu>eps % (mu>0 approx, transient case)
            errorN=max(max(abs(psin-oldpsin)))*Neps; % Use if mu is > 0.
        else % (recurrent case)
            errorN=abs(1-max(sum(psin,2)))*Neps; % Use only if mu< or = 0.
        end
        % Stopping criterion.
        if errorN<=1
            g=max_iter;
        end
        oldpsin=psin;
        g=g+1;
    end
end


end
