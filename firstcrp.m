P = [0 0;
    0 1;
    1/4 1/4;
    1 0];
v = zeros(3);

i = 0;
draw(P, i)
while i <= 0 %v(3) >= 0
    i = i + 1
    [B v] = iterate(P);
%     if v(3) >= 0
         P = B;
%     end
    draw(P, i)
end

function [P v] = iterate(P)
    % Construct A
    A = [
    P(2,1) - P(1,1), P(2,2) - P(1,2), 0 0 0 0;
    P(2,1) - P(3,1), P(2,2) - P(3,2), P(3,1) - P(2,1), P(3,2) - P(2,2), 0 0;
    0 0, P(3,1) - P(4,1), P(3,2) - P(4,2), 0 0;
    0 0, P(3,1) - P(1,1), P(3,2) - P(1,2), -1 0;
    P(2,1) - P(4,1), P(2,2) - P(4,2), 0 0 0, -1;
    ]

    % Solve (Ax = b, x = A \ b)
    B = rref(A)

    % Choose values for V
    t = 0.001;
    v = -t * B(:, 6)

    % Update P
    P(2,:) = P(2,:) + v(1 : 2)';
    P(3,:) = P(3,:) + v(3 : 4)';

    % Draw P
    P
end

function draw(P, i)
    colors = ["r" "g" "b" "c" "m" "y"];
    color = colors(1 + mod(i, length(colors)));
    hold on;
    plot([P(:, 1); P(1, 1)], [P(:, 2); P(1, 2)], "-o" + color);
    text(P(3,1), P(3,2), num2str(i));
end
