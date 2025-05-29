dx = 1E-2;
x_start = 0;
x_end = 1;

xv = x_start:dx:x_end;
nst = length(xv);

yv = zeros(1, nst);
zv = zeros(1, nst);

dydx = @(x, y, z) -2 * y + 4 * exp(-x);
dzdx = @(x, y, z) -y + z^(2/3);

xv(1) = 0;
yv(1) = 2;
zv(1) = 4;

% RK2
for i = 1 : nst - 1
    k_y1 = dydx(xv(i), yv(i), zv(i));
    k_z1 = dzdx(xv(i), yv(i), zv(i));
    k_y2 = dydx(xv(i) + dx, yv(i) + k_y1 * dx, zv(i) + k_z1 * dx);
    k_z2 = dzdx(xv(i) + dx, yv(i) + k_y1 * dx, zv(i) + k_z1 * dx);
    yv(i + 1) = yv(i) + 0.5 * (k_y1 + k_y2) * dx;
    zv(i + 1) = zv(i) + 0.5 * (k_z1 + k_z2) * dx;
end

figure(1);
hold on;
grid on;
plot(xv, yv, '-b');
plot(xv, zv, '-r');
hold off;

