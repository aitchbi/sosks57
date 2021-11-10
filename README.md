# sosks57

System of 57 spectral kernels, as show below: 

![system of 57 spectral kernels](figs/sosks57.jpg?raw=true)

which were used in e.g. the following papers for estimating the energy content of fMRI graph signals on the normalized graph Laplacian spectrum: 

> Behjat, H. and Larsson, M., 2020. Spectral characterization of functional MRI data on voxel-resolution cortical graphs. In Proc. IEEE Int.
Symp. Biomed. Imaging, 2020, pp. 558–562. [paper](https://arxiv.org/abs/1910.09507)

> Behjat, H. and Larsson, M., 2020. Characterization of spatial dynamics of fMRI data in white matter using diffusion-informed white matter harmonics. In Proc. IEEE Int.
Symp. Biomed. Imaging, 2021, pp. 1586-1590. [paper](https://doi.org/10.1101/2020.10.28.359125)

The kernels are defined on the range 0 to 2, which are the lower and upper bounds for eigenvalues of any given normalized Laplacian matrix defined as `L = I - D^{-1/2} A D^{-1/2}`, where `A` and `D` denote the graph adjacency matrix and degree matrix, respectively, and `I` denotes the identity matrix. The kernels at the lower end of the spectrum are narrower since many graph signals (in particular, fMRI voxel-wise graph signals) exhibit most of the energy content in that part of the spectrum, and thus, such narrow kernels enable more accurate estimation of signal energy profiles in this part of the spectrum. 


