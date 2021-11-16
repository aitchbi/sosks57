function [w,xtags] = spgg_get_multires_warping(sz,lmax,x,Nsbs,varargin)
% Constructs a warping function that can be used to created a systems of
% spectral kernels with narrower spectral bands at the lower end of the
% spectrum upto point x and wider bands at spectral parts above x up to
% lmax. 
%
% INPUTS: 
%   sz : 
%   lmax: 
%   x : pivot point on spectrum to transition resolution of kernels.
%   NSBs: structure, with fileds lower and total:
%       total: number of subbands to cover range [0,lmax]
%       lower: number of subbands to cover range [0,x]
%   varargin{1}: smoothen by floor(varargin{1}/sz) samples.  
%
% OUTPUTS: 
%   w : warping function, of length(xtags).   
%   xtags: spectral value of each warping value. 
%
% EXAMPLE usage: 
% Nsbs.total = 57;
% Nsbs.lower = 20;
% [w,xtags] = spgg_get_multires_warping(1e-4,2,0.1,Nsbs,0.1)
% figure;
% plot(w*lmax,xtags);
% hold on;
% d = 0:0.1:lmax;
% plot(x*ones(size(d)),d,':r','LineWidth',1)
% plot(d,(Nsbs.lower/Nsbs.total)*ones(length(d)),':r');
% set(gca,'Box','off','XLim',[0 lmax],'YLim',[0 lmax]);
%
% Hamid Behjat 
% Sep 2019.

x = x/lmax; % 0.05
y = (Nsbs.lower/Nsbs.total)*lmax; % 0.33

xtags = 0:sz:lmax; % 20001?
L = length(xtags);
w = zeros(1,L);

L1 = floor((y/lmax)*L);
L2 = L-L1-2; % works, but why 2?

d1 = 0:(x-0)/L1:x;
d2 = d1(end)+d1(2);
d3 = d2:(1-d2)/L2:1;

w(1:length(d1))     = d1;
w(length(d1)+1:end) = d3;

hf = figure;
set(hf,'Position',[10 500 500 400]); 
plot(xtags,lmax*w,'k','linewidth',2);
xlabel('original \lambda');
ylabel('warped \lambda');
grid on;

if nargin==4
    return;
end
d = floor(varargin{1}/sz);
w = smooth(w,d);
