function r = myraylrnd(rsigma, N)

if nargin < 1
    rsigma = 1;
    N = 1;
else if nargin == 1
        N = 1;
    end
end
r = sqrt(-2 * rsigma * log(myrnd(N)));

end