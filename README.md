# <span style="color:cyan">pypsc</span>: Parameter Space Concept for Determining 1D-Projected Crystal Structures

**```pypsc```** is a Python package that implements the Parameter Space Concept (PSC) for determining 1-dimensionally projected structures of independent scatterers from diffraction data. This method avoids traditional Fourier inversion and focuses on exploring all possible structural parameter combinations consistent with available diffraction data in a parameter space of dimension `m`.

## <span style="color:cyan">Features</span> 

- **No Fourier Inversion Required**: Instead of relying on Fourier sums, pypsc leverages structure factor amplitudes or intensities, represented by piece-wise analytic hyper-surfaces.
- **Isosurface Intersection**: The method obtains the coordinates of scatterers by intersecting multiple isosurfaces, allowing for the detection of all possible solutions in a single derivation.
- **Resonant Contrast**: The spatial resolution achieved may exceed traditional Fourier inversion methods, especially when resonant contrast is considered.
- **Symmetry Optimizations**: Exploits symmetry properties of isosurfaces to optimize algorithms.
- **Monte-Carlo Simulations**: Includes simulations using projections of random two- and three-atom structures to illustrate the universal applicability of the method.
- **1D Projection Efficiency**: The PSC linearization approach is more efficient than Fourier sums, working with fewer reflections.

## <span style="color:cyan">Installation</span>  

You can install `pypsc` using `pip`:

```bash
pip install git+https://github.com/mvnayagam/pypsc.git
or
pip install pypsc
```

## <span style="color:cyan">Usage</span>  

To start using `pypsc`, you can import the necessary modules and run calculations on crystal structures based on diffraction data. Though it is named `pypsc`, the modules can be called using the short form `psc`. Below is a simple example of how you might use the library:

```python
import psc
from psc.lib import hsurf_g
# so on
```
For for information, see the examples. 

## <span style="color:cyan">Methodology</span> 

### <span style="color:green">Parameter Space Concept (PSC)</span>  
The PSC method utilizes structure factor amplitudes or intensities as hyper-surfaces that define acceptable regions in parameter space. By intersecting these isosurfaces, the coordinates of scatterers are determined. This enables the detection of all possible solutions consistent with the given diffraction data.

### <span style="color:green">Key Advantages</span> 
- Detects all possible solutions from the given structure factor amplitudes.
- Potentially surpasses the spatial resolution of Fourier inversion methods.
- Exploits symmetry properties of isosurfaces for optimized computations.


## <span style="color:cyan">Examples</span>
```python
from psc.lib.g_space import g, hsurf_g
from psc.lib.x3Dlinearization import linearizenD_EPA
from psc.lib.x3Drepetition import getpolytope_EPA  
from psc.lib.x3Dchecklinearization import checklinear
from psc.lib.x3Dintersection import find_intersection
from psc.lib.x3Dreadwrite import wrtdata

# Simulate a two-atom structure and project to 1D. The parameter xcoor represents the given or assumed structure
xcoor  = np.sort(xcoor)[::-1]

# Define EPA model, where each atomic scattering factor is set to 1. 
f     = [1, 1]

# Define reflection
l = 3
# calculate amplitude for given strucutre and RO
gi    = np.abs(g(l, xcoor, f))

#---> inearization process starts here. The center of polarities are first calculated. 
meshlist = getmesh(l, xcoor, isos.max())

#---> double segment method - EPA 
pnts = double_segment_EPA(gi, l, f, error=0)    # This step does linearization of intensity gi
plist=linrep_DS(l, f, pnts, meshlist, imin=0, imax=0.5) # This step repeats the found linearized isosurface in entire PS

#---> Write the calculated linearization data in a file.
writedata(fn, plist)

# ---> Writing found solutions
analyzesolution(solution, xcoor, plotting=True)

# ----------------------------------------------------------------------------
# Plotting the actual isosurface and linearized segments.
# These routine are only available for two- and three-
# dimensional PS in present level of PSC development.
# ----------------------------------------------------------------------------

# calculate the isosurface over entire PS using above gi for s=+1 and s=-1. Ideally not necessary for pypsc. The isosurface giso1 and giso2 are just for visualization purpose.
giso1 = hsurf_g(l, grid, f, gi, j, s=1)
giso2 = hsurf_g(l, grid, f, gi, j, s=-1)

# plot calculated isosurfcae. keyword cc defines the color of isosurface
plotisosurf_EPA(l, h, gi, ax, isos, giso1, giso2, cc='C0', lw=2, imax=0.5)

#---> plot segments
plot_segment(ax, plist, cc='C0')
```

## <span style="color:cyan">Documentation</span> 

For detailed documentation on how to use the library, please visit the [Documentation](https://github.com/mvnayagam/pypsc.git).

## <span style="color:cyan">License</span> 

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## <span style="color:cyan">Contributing</span> 

Contributions are welcome! Please fork the repository, make your changes, and submit a pull request. For major changes, open an issue to discuss what you would like to change.

## <span style="color:cyan">Issues</span> 

If you encounter any problems or have suggestions, feel free to open an issue on the [Issues page](https://github.com/mvnayagam/pypsc/issues).

---