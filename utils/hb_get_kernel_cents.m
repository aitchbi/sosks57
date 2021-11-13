function y = hb_get_kernel_cents(g,sz,lmax)
% Determines centers of a given cell array of kernels.
%
%
%
% Hamid Behjat 
% Sep 2019.

y = zeros(size(g));
x = 0:sz:lmax;
for n=2:length(g)-1
    [~,d1] = max(g{n}(x));
    y(n) = x(d1);
end
y(1) = y(2)/2;
y(end) = y(end-1)+(lmax-y(end-1))/2;