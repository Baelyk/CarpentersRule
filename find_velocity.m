% Exitflags:
% 1: Using fmincon's output
% 0: Using linprog's output
% -1: Using rand output

function [velocity, exitflag] = find_velocity(polygon, is_polygon, resolution)
		% HOW TO PIN OTHER EDGES EASILY
		% Set some options:
		% - Feasibility mode helps it find a solution
		% - Increase max function evaluations and iterations because it has a hard time
		%   finding  solution
		% - Run the optimizer in parallel to use more cores
		options = optimoptions(@fmincon, ...
			"Algorithm", "sqp", ...
			"EnableFeasibilityMode", true, ...
			"Display", "none", ...
			"UseParallel", false);

		% Encode the equalities into a matrix Aeq
		Aeq = create_Aeq(polygon, is_polygon);
		beq = zeros(size(Aeq, 1), 1);
		% Encode the inequalities
		Ain = -1 * create_Ain(polygon, is_polygon);
		bin = -1 * create_bin(polygon, is_polygon, resolution);
		if options.Algorithm ~= "sqp"
			Aeq = sparse(Aeq);
			Ain = sparse(Ain);
			bin = sparse(bin);
		end
		% Set the lower and upper bound of P(1,:) and P(2,:) to 0
		fix = [0 0 0 0]';
		ub = inf(2 * size(polygon, 1), 1);
		%ub(1 : 4) = fix;
		lb = -ub;
		% Start the optimization at the zero vector, height 2*n because each point has x and y velocity
		%x0 = zeros(2 * size(polygon, 1), 1);
		% Start at a feasible solution
		f = zeros(2 * size(polygon, 1), 1);
		[x0, fval, exitflag, output, lambda] = linprog(f, Ain, bin, Aeq, beq, [], ub, optimoptions("linprog", "Algorithm", "dual-simplex", "Display", "none", "Preprocess", "basic", "Diagnostics", "off"));

		if exitflag ~= 1
			fprintf(" !! linprog failed (%d) to find a feasible point, using 0 as initial\n", exitflag);
			x0 = zeros(2 * size(polygon, 1), 1);
			%m_eq = size(Aeq, 1);
			%m_in = size(Ain, 1);
			%n = size(Aeq, 2);
			%A = sparse(m_eq + m_in, n + m_in);
			%A(1 : m_eq, 1 : n) = Aeq;
			%A(m_eq + 1 : m_eq + m_in, 1 : n) = Ain;
			%%A(m_eq + 1 : m_eq + m_in, n + 1 : n + m_in) = -eye(m_in);

			%b = [zeros(m_eq, 1); bin(:)];
			%x0 = A \ b;
			%norm(A * x0 - b)
			%x0 = full(x0(1 : n))
		end

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
			lb, ...
			... % The upper bound for the velocities of P(1) and P(2) is zero to fix them
			ub, ...
			... % donothing is a function that does nothing because we have no nonlinear
			... %   constraints, but must still provde something to the optimizer
			donothing, ...
			... % Our options, as above
			options);

		% If converged to an infeasible point, use the linear programs feasible point
		if exitflag == -2
			%if output.bestfeasible ~= []
				%fprintf(" !! fmincon converged to infeasible point, using its best feasible\n");
				%velocity = output.bestfeasible.x;
				%exitflag = 0;
			%else
				fprintf(" !! fmincon converged to infeasible point, using linprog's feasible point\n");
				velocity = x0 / norm(x0);
				%velocity(isnan(velocity)) = 0;
				exitflag = 0;
				if any(isnan(velocity))
					fprintf(" !! velocity contained NaN, using rand velocity\n");
					velocity = rand(size(velocity));
					exitflag = -1;
				end
			%end
		elseif exitflag ~= 1
			fprintf(" !! fmincon exited with %d\n", exitflag);
			exitflag = 1;
		end

	velocity = reshape(velocity, 2, [])';
end

function Aeq = create_Aeq(P, is_polygon)
	[m n] = size(P);
	Aeq = [];

	% Fix kth bar
	k = 1;
	Aeq = [zeros(4, 2 * (k - 1)) eye(4)];

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

