dx = 0.1;
x0 = 0;
xf = 5;

x = x0 : dx : xf;
nsteps = length(x);
y = zeros(1, nsteps);
u = zeros(1, nsteps);

F = @(x, y, u) -0.6 * u - 8 * y;

x(1) = 0;
y(1) = 4;
u(1) = 0;

for i = 1 : nsteps - 1
    m1 = u(i);
    k1 = F(x(i), y(i), u(i));

    m2 = u(i) + 0.5 * k1 * dx;
    k2 = F(x(i) + 0.5 * dx, y(i) + 0.5 * m1 * dx, u(i) + 0.5 * k1 * dx);

    m3 = u(i) + 0.5 * k2 * dx;
    k3 = F(x(i) + 0.5 * dx, y(i) + 0.5 * m2 * dx, u(i) + 0.5 * k2 * dx);


    m4 = u(i) + k3 * dx;
    k4 = F(x(i) * dx, y(i) + m3 * dx, u(i) + k3 * dx);

    y(i + 1) = y(i) + (m1 + 2 * m2 + 2 * m3 + m4) * (dx / 6);
    u(i + 1) = u(i) + (k1 + 2 * k2 + 2 * k3 + k4) * (dx / 6);
end

clf;
figure(1);
hold on;

title("Solution for diffeq");
xlabel("x");
ylabel("y");

plot(x, y, '-r');
plot(x, u, '-b');


legend("L $-0.6u - 8y$", "L $(-0.6u - 8y)'$", Interpreter="latex");

hold off;


