v = [
    0 3;
    2 3;
    -3 1;
    3 0;
    3 3;
    5 2;
    4 1;
    -2 -4;
    -3 -1;
    1 -2;
];
ellipse = welzl(v,[]);
c = [ellipse(1), ellipse(2)]; %Center point
a = ellipse(3); % semi-major axis
b = ellipse(4); % semi-minor axis
t = ellipse(5); % rotation
A = pi*a*b

function ellipse = welzl(interior, boundary)
    if size(interior,1) == 0
        if size(boundary,1) < 3
            ellipse = false;
            return
        elseif size(boundary,1) == 3
            ellipse = ellipse3(boundary);
            return
        elseif size(boundary,1) == 4
            ellipse = ellipse4(boundary);
            return
        else
            ellipse = ellipse5(boundary);
            return
        end
    end

    i = randi([1 size(interior,1)],1,1);
    p = interior(i,:);
    interior(i,:) = [];

    ellipse = welzl(interior,boundary);
    check = in_ellipse(p, ellipse);

    if not(check)
        boundary = [boundary; p];
        ellipse = welzl(interior,boundary);
    end
end

function check = in_ellipse(p,ellipse)
    if ellipse == false
        check = false;
        return
    end
    c = [ellipse(1), ellipse(2)];
    a = ellipse(3);
    b = ellipse(4);
    t = ellipse(5);
    p = p - c;
    rot_mat = [cos(t), sin(t); -sin(t), cos(t)];
    F = transpose(rot_mat)*diag([1/a^2,1/b^2])*rot_mat;

    if p*F*transpose(p) <= 1
        check = true;
    else
        check = false;
    end
end

function ellipse = to_geometric(F,c)
    [V,D] = eig(F);
    a = 1/sqrt(D(1,1));
    b = 1/sqrt(D(2,2));
    if V(2,1) < 0
        V(:,1) = -V(:,1);
    end
    t = acos(V(1,1));
    ellipse = [c,a,b,t];
end

function ellipse = ellipse3(boundary)
    c = mean(boundary);
    boundary = boundary - c;
    F = 1.5*inv(transpose(boundary)*boundary);
    ellipse = to_geometric(F,c);
end

function ellipse = ellipse5(S)
    x = S(:, 1);
    y = S(:, 2);
    A = [x.^2, y.^2, 2 * x .* y, x, y];

    if cond(A) >= 1 / eps
        c = [];
        a = [];
        b = [];
        t = [];
        return;
    end

    sol = linsolve(A, -ones(size(S, 1), 1));

    c = linsolve(-2 * [sol(1), sol(3); sol(3), sol(2)], sol(4:5));

    A = [eye(3), -[sol(1), sol(3), sol(3)]'; c(1)^2, 2 * c(1) * c(2), c(2)^2, -1];
    s = linsolve(A, [0; 0; 0; 1]);
    F = [s(1), s(2); s(2), s(3)];

    ellipse = to_geometric(F, c);
end

function ellipse = ellipse4(S)
    Sc = S - mean(S, 1);
    angles = atan2(Sc(:, 2), Sc(:, 1));
    [~, idx] = sort(-angles);
    S = S(idx, :);

    A = [S(3, :) - S(1, :); S(2, :) - S(4, :)];
    b = S(2, :) - S(1, :);
    s = linsolve(A, transpose(b));
    diag_intersect = S(1, :) + s(1) * (S(3, :) - S(1, :));

    S = S - diag_intersect;

    AC = S(3, :) - S(1, :);
    theta = atan2(AC(2), AC(1));
    rot_mat = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    S = S * rot_mat';

    m = (S(2, 1) - S(4, 1)) / (S(4, 2) - S(2, 2));
    shear_mat = [1, m; 0, 1];
    S = S * shear_mat';

    b = sqrt(sum(S.^2, 2));
    d = b(2) * b(4) / (b(3) * b(1));
    stretch_mat = diag([d^0.25, d^-0.25]);
    S = S * stretch_mat';

    a = sqrt(sum(S.^2, 2));
    coeff = zeros(1, 4);
    coeff(1) = -4 * a(2)^2 * a(3) * a(1);
    coeff(2) = -4 * a(2) * (a(3) - a(1)) * (a(2)^2 - a(3) * a(1));
    coeff(3) = 3 * a(2)^2 * (a(2)^2 + a(3)^2)...
        - 8 * a(2)^2 * a(3) * a(1) + 3 * (a(2)^2 + a(3)^2) * a(1)^2;
    coeff(4) = coeff(2) / 2;
    rts = roots(coeff);
    rts = rts((-1 < rts) & (rts < 1));
    theta = asin(real(rts(1)));

    D_mat = [cos(theta)^-0.5, sin(theta) * cos(theta)^-0.5; 0, cos(theta)^0.5];
    S = S * D_mat';

    boundary = S(1:3, :);
    A = [-2 * boundary, ones(3, 1)];
    b = -sum(boundary.^2, 2);
    s = linsolve(A, b);

    circle_c = s(1:2);
    circle_r = sqrt(sum(circle_c.^2) - s(3));

    T_mat = D_mat * stretch_mat * shear_mat * rot_mat;

    ellipse_c = linsolve(T_mat, circle_c)' + diag_intersect;
    ellipse_F = (T_mat' * T_mat) / circle_r^2;

    ellipse = to_geometric(ellipse_F, ellipse_c);
end