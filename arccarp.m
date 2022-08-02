% CONSTANTS
% The space between spiral arms
epsilon = 1/4;
% The distance to step with the velocity in each iteration
stepsize = 10^-3;
% E.g. (a, b +/- resolution) is considered (a, b)
resolution = stepsize;
% Estimated ticks needed, used for preallocation
ticks = 500;
% Whether P is a polygon
is_polygon = false;

% Generate the spiral
%P = spiral(epsilon)
% Keep track of all P
all_P = {};
% P is a polygon if the first and last point is the same
if P(1, :) == P(end, :)
	is_polygon = true;
	P = P(1 : end - 1, :);
end
L = spiral_length(P);
ignored = repmat(false, size(P, 1), 1);

% Clear the figure
%clf;
% Set the axis
ax = ([-epsilon, sum(L) + epsilon, -epsilon, sum(L) + epsilon]);
% Create new frames cell array
gif_filename = "spiral " + datestr(datetime, 31) + ".gif";
%draw(P, 0, ax);

% Keep track of velocities of all iterations, and preallocate
V = zeros(2 * size(P, 1), ticks);

% Set the lower and upper bound of P(1,:) and P(2,:) to 0
fix = [0 0 0 0]';
% HOW TO PIN OTHER EDGES EASILY
% Set some options:
% - Feasibility mode helps it find a solution
% - Increase max function evaluations and iterations because it has a hard time
%   finding  solution
% - Run the optimizer in parallel to use more cores
options = optimoptions(@fmincon, ...
	"EnableFeasibilityMode", true, ...
	"MaxIterations", 6e4, ...
	"MaxFunctionEvaluations", 6e4, ...
	"Display", "none", ...
	"UseParallel", false);

% Run the optimizer while the largest y-coord of any point is greater than
% epsilon. My hope was that after that it would be mostly straight.
i = 0;
q = 2;
while not_done(P, is_polygon)
	i = i + 1

    % Start the optimization at the zero vector, height 2*n because each point has x and y velocity
    %x0 = zeros(2 * size(P, 1), 1);
	% Encode the equalities into a matrix Aeq
	Aeq = create_Aeq(P, is_polygon);
	% Encode the inequalities
	Ain = -1 * create_Ain(P, is_polygon);
    bin = -1 * create_bin(P, is_polygon, resolution);
	lb = repmat(-Inf, 2 * size(P, 1), 1);
	ub = repmat(Inf, 2 * size(P, 1), 1);
	lb(end - 3 : end) = 0
	ub(end - 3 : end) = 0
	x0 = linprog(zeros(2 * size(P, 1), 1), Ain, bin, Aeq, zeros(size(Aeq, 1), 1), fix, fix, optimoptions("linprog", "Display", "iter"))

		A = eye(4);
		A(end + 1 : end + size(Aeq, 1), 1 : size(Aeq, 2)) = Aeq;
		[m, n] = size(Ain);
		A(end + 1 : end + m, 1 : n + m) = [Ain eye(m)];
		b = [zeros(4 + size(Aeq, 1), 1); bin(:)];

		B = rref([A b]);

	% This gets the optimal v.
	[v fval exitflag output] = fmincon(...
		... % @cdr_obj_fun is the objective function from the paper
		@cdr_obj_fun, ...
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
		options)
    % Store V
    V(1 : length(v), i) = v;

	v = stepsize * v;
	[P P_ignored] = update_positions(P, v, ignored, L);

    % Remove the bar stretching
    P = fix_lengths(P, L, ignored);

    % Check for colinear bars
    [P, L, ignored] = remove_colinear_bars(P, L, ignored, is_polygon, resolution);

	% Save this P
	if is_polygon
		% Add first to end if a polygon
		P_ignored(end + 1, :) = P_ignored(1, :)
	end
    all_P{i} = P_ignored;
	P

	% Draw this iteration
	% draw(P, i, ax);

	% If converged to an infeasible point, stop
	if exitflag < 0
		break
	end
end

%export_gif(gif_filename, all_P, i);
gif_filename

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
				( (P(i, :) - P(j, :)) * (v(2 * i - 1 : 2 * i) - v(2 * j - 1 : 2 * j)) ...
                    - norm(P(i, :) - P(j, :)) );
		end
	end
end

% Draw the polygon, change colors each time
function draw(P, i)
	% Cycle through the colors
    colors = ["r" "g" "b" "c" "m" "y"];
    color = colors(1 + mod(i, length(colors)));
	% Don't show the plots as they are drawn
	% Plot the arc, "-" connects the vertices, "o" puts a circle at the vertices
    plot(P(:, 1), P(:, 2), "-o" + color);
	title("i = " + num2str(i));
end

% Return empty matrices
function [c, ceq] = donothing(x)
	c = [];
	ceq = [];
end

function Aeq = create_Aeq(P, is_polygon)
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

function L = spiral_length(P)
    L = zeros(size(P, 1), 1);

    for i = 2 : size(P, 1)
        L(i - 1) = norm(P(i - 1, :) - P(i, :), 1);
    end
end

function P = fix_lengths(P, L, ignored)
    for i = 3 : length(ignored)
		if ignored(i)
			continue
		end
		offset = sum(ignored(1 : i));

        % v is the vector representing the bar from i-1 to i
        v = P(i - offset, :) - P(i - offset - 1, :);
		new_length = L(i - 1);
		j = i - 1;
		while ignored(j)
			new_length = new_length + L(j - 1);
			j = j - 1;
		end
        % w is v scaled to the length it should be without erros
        w = (new_length / norm(v)) * v;
        % delta is the difference between the two
        delta = w - v;
        % Translate all the vertices after and including i of curve by
        % delta
        P(i - offset : end, :) = P(i - offset : end, :) + delta;
    end
end

function [P, L, ignored] = remove_colinear_bars(P, L, ignored, is_polygon, resolution)
    i = 2;
    while i < size(P, 1)
        % v is the vector representing the bar from i-1 to i
        v = P(i, :) - P(i - 1, :);
        % w is the vector representing the bar from i to i+1
        w = P(i + 1, :) - P(i, :);

        if abs(-v*w' + norm(v) * norm(w)) < resolution
			% Note that we are ignoring i
			unignored = find(ignored == 0);
			ignored(unignored(i)) = true;

            % Remove the ith row of P to remove the ith vertex
            P = remove_ith_row(P, i);

            % Remove the bar stretching
            P = fix_lengths(P, L, ignored);
        else
            i = i + 1;
        end
    end

	% Only check the last-first connection if a polygon
	if is_polygon
		m = size(P, 1);
        % v is the vector representing the bar from m-1 to m
        v = P(m, :) - P(m - 1, :);
        % w is the vector representing the bar from m to 1
        w = P(1, :) - P(m, :);

        if abs(-v*w' + norm(v) * norm(w)) < resolution
			% Note that we are ignoring i
			unignored = find(ignored == 0);
			ignored(unignored(end)) = true;

            % Remove the ith row of P to remove the ith vertex
            P = remove_ith_row(P, m);

            % Remove the bar stretching
            P = fix_lengths(P, L, ignored);
        end
	end
end

function A = remove_ith_row(A, i)
    % Remove the ith row of A
    tail = A(i + 1 : end, :);
    A = A(1 : i - 1, :);
    A(i : i + size(tail, 1) - 1, :) = tail;
end

function [pos pos_all] = update_positions(pos, vel, ignored, lengths)
	% Update non-ignored positions
	vel = reshape(vel, 2, [])';
	pos = pos + vel;

	pos_all = [];
	for i = 1 : length(ignored)
		offset = sum(ignored(1 : i));
		if ignored(i)
			% If i is an "end"
			if i > max(find(ignored == 0))
				v = pos(1, :) - pos_all(i - 1, :);

				delta = (lengths(i - 1) / norm(v)) * v;
				pos_all(i, :) = pos_all(i - 1, :) + delta;
			else
				v = pos(i + 1 - offset, :) - pos_all(i - 1, :);

				delta = (lengths(i - 1) / norm(v)) * v;
				pos_all(i, :) = pos_all(i - 1, :) + delta;
			end
		else
			pos_all(i, :) = pos(i - offset, :);
		end
	end
end

function should_continue = not_done(pos, is_polygon)
	if is_polygon
		should_continue = length(convhull(pos)) - 1 ~= length(pos);
	else
		should_continue = length(pos) > 2;
	end
end
