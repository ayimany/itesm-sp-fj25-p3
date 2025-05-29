dx = 1E-5;
x_start = 0;
x_end = 2;

xv = x_start:dx:x_end;
nst = length(xv);

yv = zeros(1, nst);

dydx = @(x, y) -2 * y + x^3 * exp(-2 * x);
y = @(x) (exp(-2 * x))/(4) * (x^4 + 4);

xv(1) = 0;
yv(1) = 1;

% RK2
for i = 1 : nst - 1
    k1 = dydx(xv(i), yv(i));
    k2 = dydx(xv(i) + 0.5 * dx, yv(i) + 0.5 * k1 * dx);
    k3 = dydx(xv(i) + 0.5 * dx, yv(i) + 0.5 * k2 * dx);
    k4 = dydx(xv(i) + dx, yv(i) + k3 * dx);
    yv(i + 1) = yv(i) + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4) * dx;
end

figure(1);
hold on;
grid on;
plot(xv, yv, '-b');
hold off;

