function [t,ufe]= feuler(f, tspan, y0, Nh)
span = tspan(2)-tspan(1);
h = span/Nh;
t = tspan(1):h:tspan(2);
ufe = zeros(1,length(t));
ufe(1)=y0;
for i = 1:(length(t)-1)
    ufe(i+1)=ufe(i)+h*f(t(i),ufe(i));    
end
end