# pyPSC: Parameter Space Concept for Determining 1D-Projected Crystal Structures

**pyPSC** is a Python package that implements the Parameter Space Concept (PSC) for determining 1-dimensionally projected structures of independent scatterers from diffraction data. This method avoids traditional Fourier inversion and focuses on exploring all possible structural parameter combinations consistent with available diffraction data in a parameter space of dimension `m`.

## Features

- **No Fourier Inversion Required**: Instead of relying on Fourier sums, pyPSC leverages structure factor amplitudes or intensities, represented by piece-wise analytic hyper-surfaces.
- **Isosurface Intersection**: The method obtains the coordinates of scatterers by intersecting multiple isosurfaces, allowing for the detection of all possible solutions in a single derivation.
- **Resonant Contrast**: The spatial resolution achieved may exceed traditional Fourier inversion methods, especially when resonant contrast is considered.
- **Symmetry Optimizations**: Exploits symmetry properties of isosurfaces to optimize algorithms.
- **Monte-Carlo Simulations**: Includes simulations using projections of random two- and three-atom structures to illustrate the universal applicability of the method.
- **1D Projection Efficiency**: The PSC linearization approach is more efficient than Fourier sums, working with fewer reflections.

## Installation

You can install `pyPSC` using `pip`:

```bash
pip install git+https://github.com/mvnayagam/pyPSC.git
or
pip install pypsc
```

## Usage

To start using `pyPSC`, you can import the necessary modules and run calculations on crystal structures based on diffraction data. Below is a simple example of how you might use the library:

```python
import pypsc
```
For for information, see the examples. 

## Methodology

### Parameter Space Concept (PSC)
The PSC method utilizes structure factor amplitudes or intensities as hyper-surfaces that define acceptable regions in parameter space. By intersecting these isosurfaces, the coordinates of scatterers are determined. This enables the detection of all possible solutions consistent with the given diffraction data.

### Key Advantages:
- Detects all possible solutions from the given structure factor amplitudes.
- Potentially surpasses the spatial resolution of Fourier inversion methods.
- Exploits symmetry properties of isosurfaces for optimized computations.

### Monte-Carlo Simulations
The algorithm has been tested with Monte-Carlo simulations, illustrating its efficacy in predicting structures of random two- and three-atom models.

## Examples

### Example 04 - Two-Atom Structure 
```python
from psc.lib.g_space import g
from psc.lib.x3Dlinearizazion import linearizenD_EPA
from psc.lib.x3Drepetition import getpolytope_EPA  
from psc.lib.x3Dchecklinearization import checklinear
from psc.lib.x3Dintersection import find_intersection
from psc.lib.x3Dreadwrite import wrtdata

# Simulate a two-atom structure and project to 1D
xcoor  = np.sort(xcoor)[::-1]

# EPA model
f     = [1, 1]

# Define reflection
l = 3
# calculate amplitude for given strucutre and RO
gi    = np.abs(g(l, xcoor, f))

# calculate the isosurface over entrie PS using above gi for s=+1 and s=-1
giso1 = hsurf_g(l, grid, f, gi, j, s=1)
giso2 = hsurf_g(l, grid, f, gi, j, s=-1)

# plot calculated isosurfcae. define cc='k' if same is wanted. or isosurface colour will change automatically
r = np.random.uniform(0.0, 0.8, 3) ; cc = (r[0],r[1],r[2],1)
plotisosurf_EPA(l, h, gi, ax, isos, giso1, giso2, cc, lw=2, imax=0.5)

#---> inearization process with error of err 
errr = 0
meshlist = getmesh(l, xcoor, isos.max())

#---> double segment method - EPA
pnts = double_segment_EPA(gi, l, f, error=0)
plist=linrep_DS(l, f, pnts, meshlist, imin=0, imax=0.5)

#---> single segment method - EPA
#pnts = single_segment_EPA(gi, l, xexp, f, error=0)
#plist=linrep_SS(l, f, pnts, meshlist, imin=0, imax=0.5)

#---> plot segments
plot_segment(ax, plist, cc)

#---> storing segment data
writedata(fn, plist)

# ---> Writting found solutions from given ROs
analyzesolution(solution, xcoor, plotting=True)
```

## Documentation

For detailed documentation on how to use the library, please visit the [Documentation](https://github.com/mvnayagam/pyPSC.git).

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please fork the repository, make your changes, and submit a pull request. For major changes, open an issue to discuss what you would like to change.

## Issues

If you encounter any problems or have suggestions, feel free to open an issue on the [Issues page](https://github.com/mvnayagam/pyPSC/issues).

---

This `README.md` gives a clear overview of the project, its features, and usage. You can adapt and extend this based on the specifics of your implementation and how your library evolves. Let me know if you need more details!