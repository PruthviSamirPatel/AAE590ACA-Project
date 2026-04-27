function y = wrapTo360Local(x)
%WRAPTO360LOCAL Wrap degrees to [0, 360).
    y = mod(x, 360);
end
