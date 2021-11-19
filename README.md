# sosks57

System of 57 spectral kernels, as show below: 

![system of 57 spectral kernels](figs/sosks57.jpg?raw=true)

which were used in e.g. the following papers for estimating the energy content of fMRI graph signals on the normalized graph Laplacian spectrum: 

> Behjat, H. and Larsson, M., 2020. Spectral characterization of functional MRI data on voxel-resolution cortical graphs. In Proc. IEEE Int.
Symp. Biomed. Imaging, pp. 558–562. [paper](https://arxiv.org/abs/1910.09507)

> Behjat, H., et al., 2021. Characterization of spatial dynamics of fMRI data in white matter using diffusion-informed white matter harmonics. In Proc. IEEE Int. Symp. Biomed. Imaging, pp. 1586-1590. [paper](https://doi.org/10.1101/2020.10.28.359125)

The kernels are defined on the range 0 to 2, which are the lower and upper bounds for eigenvalues of any given normalized Laplacian matrix defined as `L = I - D^{-1/2} A D^{-1/2}`, where `A` and `D` denote the graph adjacency matrix and degree matrix, respectively, and `I` denotes the identity matrix. The kernels at the lower end of the spectrum are narrower since many graph signals (in particular, fMRI voxel-wise graph signals) exhibit most of the energy content in that part of the spectrum, and thus, such narrow kernels enable more accurate estimation of signal energy profiles in this part of the spectrum. 

The system of kernels are precomputed and made available as function handles; in particular, 57 function handles, which are saved in two files (due to the large file size) found in folder `mats` named `sosks57_lower.mat` (the first 20 kernels) and `sosks57_upper.mat` (the remaining 37 kernels). The central value (peak value) of each kernel is stored in file `sosks57_cents.mat`. 

For large graphs, graph signals can not be brought into the spectral domain via computing the Graph Fourier Transform (GFT) since it is impractical to fully diagonalize the graph Laplacian matrix. As such, it is not possible to obtain a spectral representation of a given graph signal on the graph at the resolution of eigenvalues, and therefore, the distribution of the signal energy cannot be directly computed uisng these kernels within the spectral doamin in a similar way as Frequency filtering in conventional signal processing. Nevertheless, one can *approximate* the distribution of the graph signals energy associated to each of these kernels using a polynomial approximation scheme, e.g. using Chebyshev filters as used in: 

> Hammond, D.K., et al, 2011. Wavelets on graphs via spectral graph theory. Applied and Computational Harmonic Analysis, 30(2), pp. 129-150. [paper](https://doi.org/10.1016/j.acha.2010.04.005) (see Section 6) 

To do this, one should use a suitable Chebyshev polynomial order to ensure that the approximated kernels also satisfy the Parseval frame condition. File `mats/sosks57_chebyOrds.mat` provides suitable polynomial orders for each kernel such that the resulting set of approximated system of spectral kernels do not deviate from the Parseval frame condition by no more than 0.01 at any point across the spectral range 0 to 2. 

Script `construct_multires_uniform_tight_frame.m` shows how this 57 system of spectral kernels were generated, in part using the signal-adapted tight frame design theory proposed in:

> Behjat, H., et al., 2016. Signal-adapted tight frames on graphs. IEEE Trans. Signal Process., 64(22), pp.6017-6029. [paper](https://bme.lth.se/fileadmin/biomedicalengineering/Personal_folders/Hamid_Behjat/HBehjat_TSP2016.pdf)

    
