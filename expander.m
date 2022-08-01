function [positions, velocities] = expander(polygon, stepsize)
	% E.g. (a, b +/- resolution) is considered (a, b)
	resolution = stepsize;
	% Estimated ticks needed, used for preallocation
	ticks = 500;
	% Whether P is a polygon
	is_polygon = false;

	% Keep track of all positions and velocities
	positions = {}; %cell(1, ticks);
	velocities = {}; %cell(1, ticks);
	% P is a polygon if the first and last point is the same
	if polygon(1, :) == polygon(end, :)
		is_polygon = true;
		polygon = polygon(1 : end - 1, :);
    end
	lengths = spiral_length(polygon);
	ignored = false(size(polygon, 1), 1);

	i = 0;
	while not_done(polygon, is_polygon)
		i = i + 1
		remaining_vertices = length(polygon)

		% Find the optimal solution
		[v, exitflag] = find_velocity(polygon, is_polygon, resolution);

		% Store V
		velocities{i} = v;

		% Update the polygon and position
		v = stepsize * v;
		[polygon, position] = update_positions(polygon, v, ignored, lengths);

		% Remove the bar stretching
		polygon = fix_lengths(polygon, lengths, ignored);

		% Check for colinear bars
		[polygon, lengths, ignored] = remove_colinear_bars(polygon, lengths, ignored, is_polygon, resolution);

		% Save this position
		if is_polygon
			% Add first to end if a polygon
			position(end + 1, :) = position(1, :);
		end
		positions{i} = position;

		% If converged to an infeasible point, stop
		if exitflag == -2
			break
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
	return
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

function [pos, pos_all] = update_positions(pos, vel, ignored, lengths)
	% Update non-ignored positions
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
