function y=distribution(mu,v,g,x)
x=x(:);
mu=mu(:);
v=v(:);
g=g(:);
for i=1:size(mu,1)
   d = x-mu(i);
   amp = g(i)/sqrt(2*pi*v(i));
   y(:,i) = amp*exp(-0.5 * (d.*d)/v(i));
end