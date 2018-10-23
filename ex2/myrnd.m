function [r] = myrnd(seedx, seedy, seedz)

persistent x y z;
if nargin == 3
    x = seedx;
    y = seedy;
    z = seedz;
end
if isempty(x)
    x = hex2dec('dead');
    y = hex2dec('beef');
    z = hex2dec('c0de');
end
if nargin == 1
    N = seedx;
else
    N = 1;
end


r = zeros(1, N);
for i = 1:N
    x = mod(171 * x, 30269);
    y = mod(170 * y, 30307);
    z = mod(172 * z, 30323);
    r(i) = mod(x / 30269 + y / 30307 + z / 30323, 1);
end

end