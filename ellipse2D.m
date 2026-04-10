function veb = ellipse2D(m,P,std)

[V,lambda] = eig(P);
t = linspace(0,2*pi,1000);
xs = std*sqrt(lambda(1,1))*cos(t);
ys = std*sqrt(lambda(2,2))*sin(t);
for i = 1:1000
    s = [xs(i); ys(i)];
    s = V*s;
    xs(i) = m(1) + s(1);
    ys(i) = m(2) + s(2);
end


veb = [xs; ys];


end