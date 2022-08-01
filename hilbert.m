function [pos] = hilbert(final_order)
    pos = [0.25 0.25
        0.25 0.75
        0.75 0.75
        0.75 0.25];

    for order = 2 : final_order

        pos = 0.5 * pos;
        pos = [[pos(:, 2) pos(:, 1)]; ...
            pos + [0, 0.5]; ...
            pos + [0.5, 0.5]; ...
            [0.5 - pos(:, 2) 0.5 - pos(:, 1)] + [0.5, 0]];
    end

    %clf
    %plot(pos(:, 1), pos(:, 2), "-o")
    %axis([0 1 0 1])
    %shg
end
