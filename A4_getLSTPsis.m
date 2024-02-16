% Malgorzata O'Reilly 2023.
% See the text file "instructions_and_conditions_of_use"
% for the conditions of use and how to use.

function [psin, iterationsN]=A4_getLSTPsis(Q11s,Q12s,Q21s,Q22s)

% For the stopping criterion.
Neps=1e12;
max_iter=100;

% Converting for motational convenience.
Q11=Q11s;
Q12=Q12s;
Q21=Q21s;
Q22=Q22s;
Q=[Q11,Q12;Q21,Q22];
s1=size(Q11,1);
s2=size(Q22,1);

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
        errorN=max(max(abs(psin-oldpsin)))*Neps;
        % Stopping citerion.
        if errorN<=1
            g=max_iter;
        end
        oldpsin=psin;
        g=g+1;
    end
end


end
