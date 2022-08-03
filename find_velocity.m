function [velocity, exitflag] = find_velocity(polygon, is_polygon, resolution)
		% Encode the equalities into a matrix Aeq
		Aeq = create_Aeq(polygon, is_polygon);
		beq = zeros(size(Aeq, 1), 1);
		cond(Aeq)
		% Encode the inequalities
		Ain = -1 * create_Ain(polygon, is_polygon);
		cond(Ain)
		bin = -1 * create_bin(polygon, is_polygon, resolution);
		% Set the lower and upper bound of P(1,:) and P(2,:) to 0
		fix = [0 0 0 0]';
		% Start the optimization at the zero vector, height 2*n because each point has x and y velocity
		%x0 = zeros(2 * size(polygon, 1), 1);
		% Start at a feasible solution
		size(Aeq)
		rank(Aeq)
		rank([Aeq zeros(size(Aeq, 1), 1)])
		size(Ain)
		rank(Ain)
		rank([Ain bin(:)])
		f = zeros(2 * size(polygon, 1), 1);
		[x0, fval, exitflag, output, lambda] = linprog(f, Ain, bin, Aeq, beq, fix, fix, optimoptions("linprog", "Display", "iter", "Preprocess", "basic", "Diagnostics", "on", "MaxIterations", 1e4))
		% HOW TO PIN OTHER EDGES EASILY
		% Set some options:
		% - Feasibility mode helps it find a solution
		% - Increase max function evaluations and iterations because it has a hard time
		%   finding  solution
		% - Run the optimizer in parallel to use more cores
		options = optimoptions(@fmincon, ...
			"EnableFeasibilityMode", true, ...
			"MaxFunctionEvaluations", 10e4, ...
			"MaxIterations", 10e4, ...
			"Display", "iter", ...
			"UseParallel", false);

		% This gets the optimal v.
		[velocity, ~, exitflag, output] = fmincon(...
			... % @cdr_obj_fun is the objective function from the paper
			@(v) cdr_obj_fun(v, polygon), ...
			... % x0 is our initial point, always 0
			x0, ...
			... % Ain are our inequalies
			Ain, ...
			... % Our inequalities
			bin, ...
			... % Aeq are our equalities
			Aeq, ...
			... % Our equaities are all = 0
			zeros(size(Aeq, 1), 1), ...
			... % The lower bound for the velocity of P(1) and P(2) is zero
			fix, ...
			... % The upper bound for the velocities of P(1) and P(2) is zero to fix them
			fix, ...
			... % donothing is a function that does nothing because we have no nonlinear
			... %   constraints, but must still provde something to the optimizer
			donothing, ...
			... % Our options, as above
			options);

	velocity = reshape(velocity, 2, [])';
end

function Aeq = create_Aeq(P, is_polygon)
	[m n] = size(P);
	Aeq = [];

	% Fix first and second vertices
	%Aeq = eye(4);

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
	% If this is a polygon, connect the last point to the first
	if is_polygon
		i = m;
		j = 1;
		Aeq(end + 1, 1 : 2*m) = [...
			P(j, 1) - P(i, 1), ...
			P(j, 2) - P(i, 2), ...
			zeros(1, 2 * m - 4), ...
			P(i, 1) - P(j, 1), ...
			P(i, 2) - P(j, 2), ...
		];
	end
end

function Ain = create_Ain(P, is_polygon)
	[m n] = size(P);
	Ain = [];

	% Encode the inequality between i and all vertices larger than it
	for i = 1 : m - 2
		for j = i + 2 : m
			% If polygon, the first and last are actually connected, so skip
			if i == 1 && is_polygon && j == m
				continue
			end

            x = [...
				zeros(1, 2 * (i - 1)), ...
				P(i, 1) - P(j, 1), ...
				P(i, 2) - P(j, 2), ...
                zeros(1, 2 * (j - i - 1)), ...
				P(j, 1) - P(i, 1), ...
				P(j, 2) - P(i, 2), ...
				zeros(1, 2*m - 2 * (i - 1) - 2 * (j - i - 1) - 4) ...
			];
			Ain(end + 1, 1 : 2*m) = x;
		end
	end
end

function bin = create_bin(P, is_polygon, resolution)
    m = size(P, 1);
    bin = [];

	for i = 1 : m - 2
		for j = i + 2 : m
			% If polygon, the first and last are actually connected, so skip
			if i == 1 && is_polygon && j == m
				continue
			end

			bin(end + 1) = norm(P(i, :) - P(j, :)) + resolution;
		end
	end
end

% Return empty matrices
function [c, ceq] = donothing(x)
	c = [];
	ceq = [];
end

% The objective function
function C = cdr_obj_fun(v, P)
	v;
	P;
	% m is the number of points
	m = size(P, 1);
	C = [];

	% The first sum can just be written like this
	C(end + 1) = norm(v)^2;

	% This is the second sum
	for i = 1 : m - 2
		for j = i + 2 : m
			% Any reason the i,j switch order in the norm in the paper?
			C(end + 1) = 1 / ...
				( (P(i, :) - P(j, :)) * (v(2 * i - 1 : 2 * i) - v(2 * j - 1 : 2 * j)) ...
                    - norm(P(i, :) - P(j, :)) );
		end
	end
	C;
	C = sum(C);
end

