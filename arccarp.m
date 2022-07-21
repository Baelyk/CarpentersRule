epsilon = 1/4
P = spiral(epsilon)

[m n] = size(P);
x0 = zeros(2 * size(P, 1), 1)
Aeq = sparse(create_Aeq(P))
Ain = sparse(create_Ain(P))
fix = [0 0 0 0]';
options = optimoptions(@fmincon, "EnableFeasibilityMode", true, "MaxFunctionEvaluations", 100000, "MaxIterations", 10000, "UseParallel", true);


i = 0
while max(P(:, 2)) > epsilon
	i = i + 1
	v = fmincon(@cdr_obj_fun, x0, Ain, zeros(size(Ain, 1), 1), Aeq, zeros(size(Aeq, 1), 1), [fix; realmin], fix, donothing, options)
	for j = 1 : m
		j
		P(j, :) = P(j, :) + v(2 * j - 1 : 2 * j)';
	end
	draw(P, i);
end
print("spiral " + datestr(datetime, 31), "-dpng")

function draw(P, i)
    colors = ["r" "g" "b" "c" "m" "y"];
    color = colors(1 + mod(i, length(colors)));
    hold on;
    plot(P(:, 1), P(:, 2), "-o" + color);
end

function [c, ceq] = donothing(x)
	c = []
	ceq = []
end

function Aeq = create_Aeq(P)
	[m n] = size(P);
	Aeq = [];

	% Connect each point to the next one, ignoring the first and last.
	for i = 1 : m - 1
		j = i + 1;
		Aeq(end + 1, 1 : 2*m) = [...
			zeros(1, 2 * (i - 1)), ...
			P(i, 1) - P(j, 1), ...
			P(i, 2) - P(j, 2), ...
			P(j, 1) - P(i, 1), ...
			P(j, 2) - P(i, 2), ...
			zeros(1, 2*m - 2 * (i - 1) - 4) ...
		];
	end
end

function Ain = create_Ain(P)
	[m n] = size(P);
	Ain = [];

	% Encode the inequality between i and all vertices larger than it
	for i = 1 : m - 2
		for j = i + 2 : m
			Ain(end + 1, 1 : 2*m) = [...
				zeros(1, 2 * (i - 1)), ...
				P(i, 1) - P(j, 1), ...
				P(i, 2) - P(j, 2), ...
				P(j, 1) - P(i, 1), ...
				P(j, 2) - P(i, 2), ...
				zeros(1, 2*m - 2 * (i - 1) - 4) ...
			];
		end
	end
end

function C = cdr_obj_fun(v)
	global P
	[m n] = size(P);

	C = norm(v)^2;

	for i = 1 : m - 2
		for j = i + 2 : m
			% Any reason the i,j switch order in the norm in the paper?
			C = C + 1 / ...
				((P(i, :) - P(j, :)) * (v(i : i + 1) - v(j : j + 1)) - norm(P(i, :) - P(j, :)));
		end
	end
end
