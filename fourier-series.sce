// fourier series of x^2
clc
ea = 1e-8 //absolute error
err = 1 //relative error
a = -%pi
b = %pi
n = 100
for i = 1:n
function y = f(x)
y = x^2*cos(i*x)
endfunction
c(i) = intg(a,b,f)/%pi
end
for i = 1:n
function z = g(x)
z = x^2*sin(i*x)
endfunction
d(i) = intg(a,b,g,ea,err)/%pi
end
function y1 = h(x)
y1 = x^2
endfunction
a0 = intg(a,b,h)/%pi
x = linspace(-%pi,%pi,100)
for j = 1:100
A(j) = 0
B(j) = 0
for i = 1:n
A(j) = A(j)+c(i)*cos(i*x(j))
end
for i = 1:n
B(j) = B(j)+d(i)*sin(i*x(j))
end
C(j) = a0/2+A(j)+B(j)
end
plot(x,C)
title("You are looking at the Fourier series plot of function.")
xlabel("x")
ylabel("f(x)")
