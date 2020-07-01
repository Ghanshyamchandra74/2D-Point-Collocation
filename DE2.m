%function [out] = DE([func operator order coeff])
function [de] = DE2(de)
syms x y z;
de1 = de(1,1);
de2 = de(1,2);
de3 = de(1,3);
de4 = de(1,4);
de = de4.*diff(de1,de2,de3);










