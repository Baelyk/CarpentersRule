% Generate the spiral
epsilon = 1/4;
P = spiral(epsilon)

% m is the number of points, n = 2
[m n] = size(P);

% Start the optimization at the zero vector, height 2*n because each point has x and y velocity
x0 = zeros(2 * size(P, 1), 1)
% Encode the equalities into a matrix Aeq
Aeq = sparse(create_Aeq(P))
% Encode the inequalities
Ain = sparse(create_Ain(P))
% Set the lower and upper bound of P(1,:) and P(2,:) to 0
fix = [0 0 0 0]';
% Set some options:
% - Feasibility mode helps it find a solution
% - Increase max function evaluations and iterations because it has a hard time
%   finding  solution
% - Run the optimizer in parallel to use more cores
options = optimoptions(@fmincon, ...
	"EnableFeasibilityMode", true, ...
	"MaxFunctionEvaluations", 100000, ...
	"MaxIterations", 10000, ...
	"UseParallel", true);

% Run the optimizer while the largest y-coord of any point is greater than
% epsilon. My hope was that after that it would be mostly straight.
i = 0
while max(P(:, 2)) > epsilon
	i = i + 1
	% This gets the optimal v.
	v = fmincon(...
		... % @cdr_obj_fun is the objective function from the paper
		@cdr_obj_fun, ...
		... % x0 is our initial point, always 0
		x0, ...
		... % Ain are our inequalies
		Ain, ...
		... % Our inequalities are all >= 0
		zeros(size(Ain, 1), 1), ...
		... % Aeq are our equalities
		Aeq, ...
		... % Our equaities are all = 0
		zeros(size(Aeq, 1), 1), ...
		... % The lower bound for the velocity of P(1) and P(2) is zero. Add a lower
		... %   bound for the y velocity of P(3) to the smallest positive number matlab
		... %   can represent to avoid the zero vector as a solution.
		[fix; 0; realmin], ...
		... % The upper bound for the velocities of P(1) and P(2) is zero to fix them
		fix, ...
		... % donothing is a function that does nothing because we have no nonlinear
		... %   constraints, but must still provde something to the optimizer
		donothing, ...
		... % Our options, as above
		options)

	% Update all the points based on the optimal velocities
	for j = 1 : m
		% The v(<stuff>) gets the velocity for P(j)
		P(j, :) = P(j, :) + v(2 * j - 1 : 2 * j)';
	end

	% Draw this iteration
	draw(P, i);
end
% Save the figure as a png
print("spiral " + datestr(datetime, 31), "-dpng")

% The objective function
function C = cdr_obj_fun(v)
	% This lets matlab know we're using the global P variable (our spiral)
	global P
	% m is the number of points
	[m n] = size(P);

	% The first sum can just be written like this
	C = norm(v)^2;

	% This is the second sum
	for i = 1 : m - 2
		for j = i + 2 : m
			% Any reason the i,j switch order in the norm in the paper?
			C = C + 1 / ...
				((P(i, :) - P(j, :)) * (v(i : i + 1) - v(j : j + 1)) - norm(P(i, :) - P(j, :)));
		end
	end
end

% Draw the polygon, change colors each time
function draw(P, i)
    colors = ["r" "g" "b" "c" "m" "y"];
    color = colors(1 + mod(i, length(colors)));
    hold on;
    plot(P(:, 1), P(:, 2), "-o" + color);
end

% Return empty matrices
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
