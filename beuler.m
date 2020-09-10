function [t,ube] = beuler(f, tspan, y0, Nh)

span = tspan(2)-tspan(1);
h = span/Nh;
t = tspan(1):h:tspan(2);
ube = zeros(1,length(t));
ube(1)=y0;

for i = 1:(length(t)-1)
    ube(i+1)=ube(i)+h*f(t(i+1),ube(i+1));
end
end