% Part 2

mu = 1/82.45;
r0 = [-mu; 0];
r1 = [1 - mu; 0];
B = [0, 1;-1, 0];

drrdt= @(r, dr) -(1-mu) * (r - r0) / norm(r - r1)^3 + 2 * B * dr;
drdt = @(rx, ry, drx, dry) [
    drx;
    dry;
    -(1-mu) * (rx - r0(1))
];