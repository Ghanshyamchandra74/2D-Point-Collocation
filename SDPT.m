function [eqn] = SDPT(K)
%2D Langrangian Approximating Function Generating EQN
% There are two approaches Langrangian and serendipity approach 
% from which serendipity is used for less computation time 
% while the langrangian is Used for its higher Accuracy 
% Hence This Function can be used extensively where No of 
% Elements are less or Single Element Such as Weighted REsidual 
% Techniques 
% Ghanshyam_Chandra_ME_NITRR

K1=K*4; %K1 = Sum of No of Elements
% SYMBOLIC CONSTANTS INITIATION 
A = sym('c',[K1,1]);
syms x y;eqn = A(1,1)+A(2,1)*x+A(3,1)*y+A(4,1)*(x)*(y);j=4;
for k=2:K
    j=j+1;
eqn = eqn +A(j,1)*(x^k)*y+A(j+1,1)*x*(y^k)+A(j+2,1)*(x^k)+A(j+3,1)*(y^k);j=j+3;
end
end

