function [eqn] = LGN(K)
%2D Langrangian Approximating Function Generating EQN
% There are two approaches Langrangian and serendipity approach 
% from which serendipity is used for less computation time 
% while the langrangian is Used for its higher Accuracy 
% Hence This Function can be used extensively where No of 
% Elements are less or Single Element Such as Weighted REsidual 
% Techniques 
% Ghanshyam_Chandra_ME_NITRR
K1=1; %K1 = Sum of No of Elements
for i=1:K
    K1=K1+(2*i+1);
end
% SYMBOLIC CONSTANTS INITIATION 
A = sym('c',[K1,1]);
syms x y;eqn = A(1,1);j=1;
for k=1:K
for i=1:k
    j=j+1;
eqn = eqn +A(j,1)*(x^k)*(y^(k-i))+A(j+1,1)*(x^(k-i))*(y^k);j=j+1;
end
eqn = eqn + A(j+1,1)*(x^k)*(y^k);j=j+1;
end

                                                                                                                                                                                                                      
