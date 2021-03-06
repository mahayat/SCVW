function [t,ucn] = cranknic(f, tspan, y0, Nh)

span = tspan(2)-tspan(1);
h = span/Nh;
t = tspan(1):h:tspan(2);
ucn = zeros(1,length(t));
ucn(1)=y0;

for i = 1:(length(t)-1)
     new_f=@(x) ucn(i)+0.5*h*(f(t(i),ucn(i))+f(t(i+1),x))-x;
     ucn(i+1)= fzero(new_f,ucn(i));
end
end