sz         = 1e-4;
lmax       = 2;
sOrder     = 3;
Nsbs.total = 57;
Nsbs.lower = 20;
xPivot     = 0.1; 
Ls         = 0.1; % spectral length of smoothing
maxOrd     = 800; % maximum polynomial order
minOrd     = 10;  % minimum polynomial order

arange = [0 lmax];
xxx    = 0:sz:lmax;

addpath(genpath('utils')); 

%-1.Construct suitable warping
[warping,xtags] = spgg_get_multires_warping(sz,lmax,xPivot,Nsbs,Ls);

%-2.Construct signal-adapted system of spectral kernels
% [takes a minute or so]
[g,cents] = spgg_filter_design(lmax,Nsbs.total,...
    'designtype','signal_adapted_spline_type',...
    'pou','over2ndPower',...
    'sOrder',sOrder,...
    'sz',sz,...
    'warping',warping,...
    'E',xxx);

%-3.Estimate required Chebyshev polynomial orders to approximate each kernel
% [non-parallel: around 5 mins]
% [parallel, 8-Core Intel Core i9: around 90 secs]
opts=struct;
opts.maxOrds = ones(size(g))*maxOrd;
opts.minOrds = ones(size(g))*minOrd;
opts.tol.kernel     = 1e-2; % max approximation error tolerance for each kernel
opts.tol.tightframe = 1e-2; % max tight frame deviation tolerance
opts.sz = sz;
opts.parallelize = false;
[cOrds,e] = spgg_cheby_order_est(g,arange,opts);

%-4.Plot system of kernels. 
hf = figure;
set(hf,'position',[500 1000 1500 500]);
subplot(211)
spgg_view_design(g,[0,lmax],...
    'eigsToInterp',xxx,...
    'chebyOrder',cOrds.tightframe,...
    'kernelCenters',cents,...
    'tightFrameBounds',tol.tightframe);
title('System of 57 spline-type spectral kernels')
set(gca,'Box','off','XLim',[0 lmax],'YLim',[0 1.1]);
subplot(212)
plot(warping*lmax,xtags);
hold on;
d = 0:0.1:lmax;
plot(0.1*ones(size(d)),d,':r','LineWidth',1)
title('Warping function used to split spectrum into a lower & upper part')
set(gca,'Box','off','XLim',[0 lmax],'YLim',[0 lmax]);
