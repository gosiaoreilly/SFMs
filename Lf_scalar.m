% Malgorzata O'Reilly 2023.
% See the text file "instructions_and_conditions_of_use"
% for the conditions of use and how to use.

% We compute the probability density psi(t)_{ij} that
% the SFM first returns to level 0 at time t and does in phase j,
% given start from level 0 in phase i.

function FNs_mat=Lf_scalar(s_mat,phasei,phasej)

load examplepar.mat


for u=1:size(s_mat,1)
    for v=1:size(s_mat,2)
        ss=s_mat(u,v);
        X=real(ss);
        Y=imag(ss);
        
        % This is for an example with nonzero rates c_i.
        if s0==0
            Q11s=T11-ss*eye(s1);
            Q22s=T22-ss*eye(s2);
            Q12s=T12;
            Q21s=T21;
        end
        % This is for an example with zero rates c_i.
        if s0>0
            Q11s=inv(C1)*(T11-ss*eye(s1)-T10*inv(T00-ss*eye(s0))*T01);
            Q22s=inv(-C2)*(T22-ss*eye(s2)-T20*inv(T00-ss*eye(s0))*T02);
            Q12s=inv(C1)*(T12-T10*inv(T00-ss*eye(s0))*T02);
            Q21s=inv(-C2)*(T21-T20*inv(T00-ss*eye(s0))*T01);
        end
        
        [psin, iterationsN]=A4_getLSTPsis(Q11s,Q12s,Q21s,Q22s);
        
        [FNs]=psin(phasei,phasej);
        
        FNs_mat(u,v)=FNs;
        
    end
end

end
