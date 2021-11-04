clc,clf
e = 1.6e-19
k=9e9
a0 = 5e-11
smin = 0.001
smax = 10
l = input('Value for l: ')
n = input('Grid interval value: ');
s = linspace(smin,smax,n);
h = s(2) - s(1)
//potential matrix generation
V = zeros(n, n)
for j = 1:n
 V(j,j) = 2/s(j) - l*(l+1)/(s(j)*s(j))
end
//K.E. matrix generation
K = eye(n, n)*(-2)
for j = 1:n-1
 K(j,j+1)=1
 K(j+1,j)=1
end
//hamiltonian matrix generation
H=-((k*e*e)/(2*a0*h*h))*k-((k*e*e)/(2*a0))*v;
[U,EV] = spec(H);
E = diag(EV)/e;
disp("Ground state energy: "+string(E(2))+"eV");
disp("1st excited state energy: "+string(E(3))+"eV");
//experimental wavefunction
subplot(2,1,1)
plot(s',[U(:,2)],'linwidth',3);
xlabel("s(r/a0)")
ylabel("wavefunction")
legend("experimental")
//theoretical wavefunction
subplot(2,1,2)

function output=g(input)
 output = (%pi^-0.5)*(a0^-1.5)*exp(-s)
endfunction
plot(s,g(s));
xlabel("s(r/a0)")
ylabel("wavefunction")
legend("Theoretical")
