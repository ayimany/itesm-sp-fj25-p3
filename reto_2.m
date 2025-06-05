clear all, close all
m0 = 4*pi*1e-7;
x = linspace(-1, 1, 20);
y = linspace(-1, 1, 20);
z = linspace(-1, 1, 20);
[X,Y,Z] = meshgrid(x,y,z);

side_length = 0.8;
I = 5;
wires = [
    -side_length/2, -side_length/2, 0,  side_length/2, -side_length/2, 0;
     side_length/2, -side_length/2, 0,  side_length/2,  side_length/2, 0;
     side_length/2,  side_length/2, 0, -side_length/2,  side_length/2, 0;
    -side_length/2,  side_length/2, 0, -side_length/2, -side_length/2, 0
];

Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));

for i = 1:numel(X)
    B_total = [0, 0, 0];
    for w = 1:size(wires, 1)
        lvinit = wires(w, 1:3);
        lvend = wires(w, 4:6);
        a = lvinit - lvend;
        b = lvinit - [X(i), Y(i), Z(i)];
        c = lvend - [X(i), Y(i), Z(i)];
        cca = cross(a,b);
        if norm(cca) < 1e-10
            B = [0, 0, 0];
        else
            B = (m0*I/(4*pi)) * (cca/norm(cca)^2) * (dot(a,c)/norm(c) - dot(a,b)/norm(b));
        end
        B(isnan(B)) = 0;
        B(isinf(B)) = 0;
        B_total = B_total + B;
    end
    Bx(i) = B_total(1);
    By(i) = B_total(2);
    Bz(i) = B_total(3);
end

quiver3(X,Y,Z,Bx,By,Bz)

m0 = 4*pi*1e-7;
x = linspace(-1, 1, 20);
y = linspace(-1, 1, 20);
z = -0.5;
[X,Y,Z] = meshgrid(x,y,z);

side_length = 0.8;
I = 5;
wires = [
    -side_length/2, -side_length/2, 0,  side_length/2, -side_length/2, 0;
     side_length/2, -side_length/2, 0,  side_length/2,  side_length/2, 0;
     side_length/2,  side_length/2, 0, -side_length/2,  side_length/2, 0;
    -side_length/2,  side_length/2, 0, -side_length/2, -side_length/2, 0
];

Bx = zeros(size(X));
By = zeros(size(Y));
Bz = zeros(size(Z));

for i = 1:numel(X)
    B_total = [0, 0, 0];
    for w = 1:size(wires, 1)
        lvinit = wires(w, 1:3);
        lvend = wires(w, 4:6);
        a = lvinit - lvend;
        b = lvinit - [X(i), Y(i), Z(i)];
        c = lvend - [X(i), Y(i), Z(i)];
        cca = cross(a,b);
        if norm(cca) < 1e-10
            B = [0, 0, 0];
        else
            B = (m0*I/(4*pi)) * (cca/norm(cca)^2) * (dot(a,c)/norm(c) - dot(a,b)/norm(b));
        end
        B(isnan(B)) = 0;
        B(isinf(B)) = 0;
        B_total = B_total + B;
    end
    Bx(i) = B_total(1);
    By(i) = B_total(2);
    Bz(i) = B_total(3);
end

quiver3(X,Y,Z,Bx,By,Bz)

