function r = mygaussrnd(rsigma, N)

if nargin < 1
    rsigma = 1;
    N = 1;
else if nargin == 1
        N = 1;
    end
end
xr = myraylrnd(rsigma, N);
r = xr .* cos(2 * pi * myrnd(N));

end