# 2D4::xautocorr
This is an advanced pattern finder for 2D dislocation configurations. Such configurations are created by 2D4::init_config_gen and used by 2D4:.2D_DDD_simulation.

## Dependencies
* FFTW
* boost/program_options

## Parameters
* `--input-path` or `-i`: **mandatory argument**, the input filename containing the dislocation configuration (ending in `.dconf`) or the file containging the input filenames (ending in `.ini`) 1 filename per line.
* `--resolution` or `-r` (1024): the Fourier transform and autocorrelation will be calculated in this many points. (The Fourier will contain 1 extra value.) 
* `--method` or `-m`: how the dislocation density should be calulated. There are 3 different possibilities:
  * `bc`: box counting. A mesh with resolution × resolution number of cells are created and in each cell the dislocation densities are calculated based on how many dislocations can be found inside the cell.
  * `gc`: Gauss convolution. Each dislocation is represented with a Gauss-distribution of half-width `-half-width`. The distribution is a periodic repetition of a Gauss-distributions instead of sum of periodic images, therefore, it is not an exact peroidic convolution (e.g. the derivatives are not continuous), nor normalized. This is not a problem if half-width ≪ 1. The densities are sampled in the centerpoint of the mesh and due to the discrete sampling, the distribution is again not normalized, it is ensured after the calculation: ρ₊ and ρ₋ are calculated separatedly and normalized, then comes the calculation of ρ<sub>t</sub> and κ.
  * `wspn`: Wigner-Seitz positive and negative. The simulation space is meshed up to domains based on the Wigner-Seitz cell structure of the set of positive and negative dislocations. The positive ρ₊ and negative ρ₋ dislication density is then 1 / the size of the cell. The total dislocation density ρ<sub>t</sub> and signed dislocation κ density are the derived quantities, calculated by the sum and difference of these maps.
  * `wsts`: Wigner-Seitz total and signed. Like `wspn`, but the total dislocation density ρ<sub>t</sub> and signed dislocation density κ will be calculated directly and ρ₊ and ρ₋ are the derived quantities.
* `--sub-sampling` or `-s` (1): increases the size of the mesh to calculate the densities, but the resolution remains the same. Instead of calculating the densities on a mesh with size resolution × resolution, it is calculated on a mesh with size (s·r) × (s·r). This paramterer takes effect only for Gauss convolution and Wigner-Seitz methods.
* `--half-width` or `-w` (0.125): half-width of the Gauss distribution used in the Gauss convolution method. The hald-width is measured in the units of the simulation space units.
* `--output-foldername` or `-o` ("xautocorr"): the output foldername relative to the executable.

## Methods
### Box counting method

## Plot the results
### Density maps
