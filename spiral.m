function P = spiral(epsilon)
	P = [0 0;
		1 0;
		1 1;
		0 1;
		0 epsilon];

	while true
        if P(end - 2, 1) - P(end - 1, 1) < 2 * epsilon
			break
		end
		P(end + 1, :) = P(end - 3, :) + [-1, 1] * epsilon;

		if P(end - 2, 2) - P(end - 1, 2) < 2 * epsilon
			break
		end
		P(end + 1, :) = P(end - 3, :) + [-1, -1] * epsilon;

		if P(end, 1) - P(end - 2, 1) < 2 * epsilon
			break
		end
		P(end + 1, :) = P(end - 3, :) + [1, -1] * epsilon;

		if P(end, 2) - P(end - 2, 2) < 2 * epsilon
			break
		end
		P(end + 1, :) = P(end - 3, :) + [1, 1] * epsilon;
    end
end