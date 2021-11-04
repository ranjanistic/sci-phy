clc;
clf
a0=8e-16;
h=1.054e-34;
w=9.78e22;
m=1.67e-27;
smin=1e-6;
smax=5;
l=input("Enter the angular orbital number l: ")
n=input("Enter the grid interval n: ")
c1=0;
c2=10;
c3=30;
b1=c1*(1.6e32)
b2=c2*(1.6e32)
b3=c3*(1.6e32)
r1=(2*m*b1*((a0)^5))/(3*(h^2))
r2=(2*m*b2*((a0)^5))/(3*(h^2))
r3=(2*m*b3*((a0)^5))/(3*(h^2))
s=linspace(smin,smax,n)
d=s(2)-s(1)
//Creation of Potential Energy Matrix//
V1=zeros(n,n)
for j=1:n
 V1(j,j)=((s(j))^2)+(r1*((s(j))^3))+((l*(l+1))/(s(j)^2))
end
V2=zeros(n,n)
for j=1:n
 V2(j,j)=((s(j))^2)+(r2*((s(j))^3))+((l*(l+1))/(s(j)^2))
end
V3=zeros(n,n)
for j=1:n
 V3(j,j)=((s(j))^2)+(r3*((s(j))^3))+((l*(l+1))/(s(j)^2))
end
//Creation of Kinetic Energy Matrix//
K=eye(n,n)*(-2)
for j=1:n-1
 K(j,j+1)=1
 K(j+1,j)=1
end
//Creation of Hamiltonian Matrix//
H1=-(((h*w)/(2*d*d))*K)+(((h*w)/2)*V1)
H2=-(((h*w)/(2*d*d))*K)+(((h*w)/2)*V2)
H3=-(((h*w)/(2*d*d))*K)+(((h*w)/2)*V3)
[U1,EV1]=spec(H1)
[U2,EV2]=spec(H2)
[U3,EV3]=spec(H3)
E1=(diag(EV1))/((1e6)*(1.6e-19))
disp("Ground state energy for b1=0 is "+string(E1(1))+"MeV")
E2=(diag(EV2))/((1e6)*(1.6e-19))
disp("Ground state energy for b2=10 is "+string(E2(1))+"MeV")
E3=(diag(EV3))/((1e6)*(1.6e-19))
disp("Ground state energy for b3= 30 is "+string(E3(1))+"MeV")
//Experimental Wavefunction Plot//
title("Ground state wave function of anharmonic oscillator")
subplot(3,1,1)
plot(s',[U1(:,1)],'red');
xlabel("s(r/a0)")
ylabel("Wavefunction")
legend('Ground State for b1=0')
subplot(3,1,2)
plot(s',[U2(:,1)],'blue');
xlabel("s(r/a0)")
ylabel("Wavefunction")
legend('Ground State for b2=10')
subplot(3,1,3)
plot(s',[U3(:,1)],'black');
xlabel("s(r/a0)")
ylabel("Wavefunction")
legend('Ground State for b3=30')
