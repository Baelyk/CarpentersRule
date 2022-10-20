function P = spiral(epsilon)
	P = [0 0;
		1 0;
		1 1;
		0 1;
		0 epsilon];

	points = 2 + 2 / epsilon - 5;

	for i = 1 : points
		point = mod(i, 4);
        if point == 1
			P(end + 1, :) = P(end - 3, :) + [-1, 1] * epsilon;
		end

        if point == 2
			P(end + 1, :) = P(end - 3, :) + [-1, -1] * epsilon;
		end

        if point == 3
			P(end + 1, :) = P(end - 3, :) + [1, -1] * epsilon;
		end

        if point == 0
			P(end + 1, :) = P(end - 3, :) + [1, 1] * epsilon;
		end
    end
end
