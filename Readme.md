| arXiv |
|:-----:|
|[![arXiv](https://img.shields.io/badge/arXiv-2511.04735-orange.svg)](https://arXiv.org/abs/2511.04735)|

# **NuFast-Earth**: A fast code for neutrino oscillations through the Earth from the atmosphere, the Sun, or a supernova

### Overview
**NuFast-Earth** is designed to calculate all nine neutrino oscillation probabilities in matter through the Earth for atmospheric, solar, and supernova neutrinos sources. The code is designed to be very fast and it also automatically reuses repeated calculations when possible.

### General usage
For details on the usage, please read the included `Guide.pdf` file as well as the associated paper on the arXiv.

Not all files are necessary for use. In particular, the `examples` subdirectories containing `Figures.cpp`, `Speed.cpp`, and `Validation.cpp` and their associated header files are not necessary.

If one is interested in long-baseline calculations through approximately constant density profiles, we encourage the use of the faster `NuFast-LBL` code available at [https://github.com/PeterDenton/NuFast-LBL](https://github.com/PeterDenton/NuFast-LBL).

We strongly encourage a careful examination of the constants used (G<sub>F</sub>, etc.) when comparing to existing codes.

### Usage
If you use this code, please cite the associated paper [arXiv:2511.04735](https://arxiv.org/abs/2511.04735) by Peter Denton and Stephen Parke. Please also let us know if you find any bugs or further optimizations or if you run your own speed tests.
