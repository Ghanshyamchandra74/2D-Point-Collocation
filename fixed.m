%2D_point_collocation

tic;clear;clc;KAF = 5;ns=2;[e] = LGN(KAF);pretty(e);n=3;R=0;syms x y;
%***********************************************
% Value of n & m for profiles 

Lx=60;Lxx=Lx;Ly=1;Lyy=Ly;a=100/2;



w = 25;%25000;%10 25 50
KK = 0.25;%;*(((Lxx-a)/Lyy)^2);% 1 0.5 0.25
%*********************************************

hs =10.45+10*(x*w*Lyy*(1e-03))^(0.5)-(x*w*Lyy*(1e-03));
Kr=200;

K1=1; %K1 = Sum of No of Elements
for i=1:KAF
    K1=K1+(2*i+1);
end
f1=K1;

Bis = (hs/(((KK)^2)*Kr))*(Lyy*(1e-03));
Bic = (1);%(0.1)*((Lxx-a)/Lyy);
%CRF = input('Convective Rotation Factor = ');w = input('Angular Velocity of Disc = ');h0 = input('Satatic Convevctive Heat Transfer Coefficient =  ');
%hs = (CRF*x*w+h0);
r=0;
    %fprintf('For %f deqn\n',i);
    A(1,1) = x; A(2,1) = x; A(3,1) = y; 
    
    
    B(1,1)=1;B(2,1)=2;B(3,1)=2;
    C(1,1)=(1/x);C(2,1)=1;C(3,1)=KK*KK;
    for i=1:n
    R = R+(DE2([e A(i,1) B(i,1) C(i,1)]));
end 
R=R+r;EQN=e+r;b=2;

    %fprintf('For %f deqn\n',i);
    D(1,1) = x;D(2,1) = x;
    E(1,1)=1;E(2,1)=1;
    F(1,1)=1;F(2,1)=1;
    r1 = -Bic*(1-e);
    %for i=1:b
    G(1,1) = (DE2([EQN D(1,1) E(1,1) F(1,1)]))-r1;G(2,1) = (DE2([EQN D(2,1) E(2,1) F(2,1)]))-0;
%end
%disp('no of Boundary conditions for Y - direction  = ');
%for i=1:b
    %fprintf('For %d deqn\n',i);
    D(3,1) = y;  D(4,1) = y; 
    E(3,1)=1; E(4,1)=1;
    F(3,1)=1;F(4,1)=1;
    r2 = -Bis*e;
    G(3,1) = (DE2([EQN D(3,1) E(3,1) F(3,1)]))-0;G(4,1) = (DE2([EQN D(4,1) E(4,1) F(4,1)]))-r2;
%end
%if Lx>=Ly 
    L = Lx;point=ns;
    Lx = linspace(a,Lxx,point+2);Lx=Lx';Ly = linspace(0,Lyy,point+2);Ly=Ly';
    for i=1:b
    if mod(i,2)==1
    for j=1:point
    H(j,1) = subs(G(i,1),{x y},{a Ly(j,1)});
    end
    for j = 1:point
        K(j,1) = subs(G(i+b,1),{x y},{Lx(j+1,1) 0});
    end
    else
        for j=1:point
    H(j+point,1) = subs(G(i,1),{x y},{Lxx Ly(j+1,1)});
    end
    for j = 1:point
        K(j+point,1) = subs(G(i+b,1),{x y},{Lx(j+1,1) Lyy});
    end
  end
end
%else
   % L=Lyy;point=input('No of points to be satisfied in Y - Direction = ');
    %Ly = linspace(0,Lyy,point+2);Ly=Ly';Lx = linspace(a,Lxx,point+2);Lx=Lx';
    %for i=1:b
    %if mod(i,2)==1
    %for j=1:point
    %H(j,1) = subs(G(i,1),{x y},{a Ly(j+1,1)});
    %end
    %for j = 1:point
     %   K(j,1) = subs(G(i+b,1),{x y},{Lx(j+1,1) 0});
    %end
    %else
     %   for j=1:point
    %H(j+point,1) = subs(G(i,1),{x y},{Lxx Ly(j+1,1)});
    %end
    %for j = 1:point
     %   K(j+point-1,1) = subs(G(i+b,1),{x y},{Lx(j+1,1) Lyy});
    %end
  %end
%end
%end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




f2=f1-(2*point+2*point);
if mod(f2,4) == 0
f2x=linspace(a,Lxx,point+f2);f2x=f2x';f2y=linspace(0,Lyy,point+f2);f2y=f2y';
    for i=1:f2/4 
    J(i,1) = subs(R,{x y},{a f2y(i+1,1)});
    end

   i11 = length(J);
    
    for i=1:f2/4
    J(i+i11,1) = subs(R,{x y},{Lxx f2y(i+1,1)});
    end
      i12 = length(J);
    
    for i=1:f2/4
    J(i+i12,1) = subs(R,{x y},{f2x(i+1,1) 0});
    end
      i13 = length(J);
    for i=1:f2/4
    J(i13+i,1) = subs(R,{x y},{f2x(i+1,1) Lyy});
    end
else
    f2x=linspace(a,Lxx,f2*6);f2x=f2x';f2y=linspace(0,Lyy,f2*6);f2y=f2y';
    for i=1:f2
    J(i,1) = subs(R,{x y},{f2x(i+1,1) f2y(i+1,1)});
    end
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if J==0
 [Z V] = equationsToMatrix([H; K;],[sym('c',[1,f1])]);
else
[Z V] = equationsToMatrix([H; K; J;],[sym('c',[1,f1])]);
end
     disp('Matrix of coeffs = '); disp(Z);disp('Value of constants ='); disp(V);disp('Not a Unique Solution ');
     P = linsolve(Z,V);
disp('Matrix of coeffs = ');disp(Z);disp('Value of constants =');disp(V);
syms F %eqn formation for final field variable %Initialize f(x)=0 %Formation of Field Variable
K1=1; %K1 = Sum of No of Elements
for i=1:KAF
    K1=K1+(2*i+1);
end
% SYMBOLIC CONSTANTS INITIATION 
A = sym('c',[K1,1]);
syms x y;F = P(1,1);j=1;
for k=1:KAF
for i=1:k
    j=j+1;
F = F +P(j,1)*(x^k)*(y^(k-i))+P(j+1,1)*(x^(k-i))*(y^k);j=j+1;
end
F = F+ P(j+1,1)*(x^k)*(y^k);j=j+1;
end
fprintf('Field Variabele = ');
disp(F);
poin= 50;   %input('soln at how many points(excluding boundary) = ');
sd=   15; %input('steps in distribution of field variable =')
t1=linspace(a,Lxx,poin+2);t2=linspace(0,Lyy,poin+2);t1=t1';t2=t2';t21=linspace(0,-Lyy,poin+2);t21=t21';t11=linspace(a,Lxx,poin+2);t22=linspace(Lyy,0,poin+2);


for i=1:poin+2
C1(i,1) = subs(F,{x y},{t1(i,1) t2(i,1)});
C2(i,1) = subs(F,{x y},{t1(i,1) t21(i,1)});
%fprintf('field variable at {x y} = {%f %f} ===> %f \n',t1(i,1),t2(i,1),C1(i,1));
end
[xx yy] = meshgrid(t1,t2);[xx yy1] = meshgrid(t1,t21);[xx2 yy2] = meshgrid(t11,t22);zz = subs(F,{x y},{xx2 yy});clf;

if Lxx>=Lyy
figure(1);[C3,h1] = contourf(xx,yy,zz,sd);clabel(C3,h1); hold on;
else
figure(1);[C3,h1] = contourf(yy,xx,zz,sd);clabel(C3,h1); hold on;    
end

figure(1);[C2,h2] = contourf(xx,yy1,zz,sd); hold on;colorbar;colormap jet;
title('Temperature Distribution (Normalised)');
xlabel('Radius(Normalised)');
ylabel('Thickness(mm)');
sd2=50;
sd3=10;
[r,theta]=meshgrid(t1,(linspace(0,2*pi,length(t1)))');
Temp = subs(F,{x y},{r meshgrid(t2)});
figure(2);contourf(r.*cos(theta),r.*sin(theta),Temp,sd2);colorbar;colormap jet;
title('Temperature Distribution {Axis Symmetric}');
xlabel('Radius(m)');
ylabel('Thickness(m)');
[r,theta]=meshgrid(t1,(linspace(0,pi,length(t1)))');
figure(3);[c5 h5]=contourf(r.*cos(theta),r.*sin(theta),Temp,sd3);clabel(c5,h5);colorbar;colormap jet;
title('Temperature Distribution {Top View}');
xlabel('Radius(m)');
ylabel('Thickness(m)');

%*********************************

%hs=subs(hs,{x w},{a w});
%Calculation Of Heat 
%disp('At Suraface using Convective');
%fprintf('%f \n',( 2*hs*(pi)*((Lxx^2-a^2)*50*subs(F,{x y},{a Lyy})))*(1e-06));
%fprintf('%f',( 4*pi*2000*a*Lyy*50*(1-subs(F,{x y},{0 0}))));
%fprintf('%f \n',-( 4*pi*a*KK*50*subs(diff(F,x,1),{x y},{a 0})));
%******************************
%toc;

%fprintf('%f \n',(subs(F,{x y},{t1(:,1) t2(:,1)*0})))



%Efficiency of fin
%Nfin = -((Lyy*a*(subs((diff(F,x)),{x y},{a 0})))/((KK^2)*Bis*((Lxx)^2-a^2)));
%fprintf('%f \n',subs(Nfin,t222,t2));
%fprintf('%f \n',subs(Nfin,x,t1));
%Effectiveness of fin
%disp('Effectiveness\n');
%Eff = Nfin*(((Lxx)^2-a^2))/(a);
%fprintf('%f \n',subs(Eff,x,t1));
%plot(t1,subs(Nfin,x,t1));