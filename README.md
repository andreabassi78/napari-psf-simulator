# napari-psf-simulator

[![License](https://img.shields.io/pypi/l/napari-psf-simulator.svg?color=green)](https://github.com/andreabassi78/napari-psf-simulator/raw/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/napari-psf-simulator.svg?color=green)](https://pypi.org/project/napari-psf-simulator)
[![Python Version](https://img.shields.io/pypi/pyversions/napari-psf-simulator.svg?color=green)](https://python.org)
[![tests](https://github.com/andreabassi78/napari-psf-simulator/workflows/tests/badge.svg)](https://github.com/andreabassi78/napari-psf-simulator/actions)
[![codecov](https://codecov.io/gh/andreabassi78/napari-psf-simulator/branch/main/graph/badge.svg)](https://codecov.io/gh/andreabassi78/napari-psf-simulator)
[![napari hub](https://img.shields.io/endpoint?url=https://api.napari-hub.org/shields/napari-psf-simulator)](https://napari-hub.org/plugins/napari-psf-simulator)

A plugin for the simulation of the 3D Point Spread Function of an optical systen, particularly a microscope objective.
 
Calculates the PSF using scalar and vectorial models.  
The following aberrations are included:
- phase aberration described by a Zernike polynomials with n-m coefficients.
- aberration induced by a slab, with a refractive index different from the one at the object.  

----------------------------------

This [napari] plugin was generated with [Cookiecutter] using [@napari]'s [cookiecutter-napari-plugin] template.

<!--
Don't miss the full getting started guide to set up your new package:
https://github.com/napari/cookiecutter-napari-plugin#getting-started

and review the napari docs for plugin developers:
https://napari.org/plugins/stable/index.html
-->

## Installation

You can install `napari-psf-simulator` via [pip]:

    pip install napari-psf-simulator


To install latest development version :

    pip install git+https://github.com/andreabassi78/napari-psf-simulator.git


## Usage

1) Lauch the plugin and select the parameters of the microscope: `NA` (numerical aperture), `wavelenght`, `n` (refractive index at the object),
   `FOV xy` (field of view in the transverse direction), `FOV z` (field of view in the axial direction), `dxy` (pixel size, transverse sampling), `dz` (voxel depth, axial sampling), `lens radius` (physical aperture of the lens, used in vectorial model)

2) Select a propagation model between `scalar` and `vectorial`.  

3) Select an aberration type (if needed).

4) Press `Calculate PSF` to run the simulator. This will create a new image layer with the 3D PSF.
 
   The option `Show Airy disk` creates a circle with radius given by the diffraction limit (Rayleigh criterion).

![raw](https://github.com/andreabassi78/napari-psf-simulator/raw/main/images/figure.png)
**Napari viewer with the psf-simulator widget showing the in-focus plane of an aberrated PSF**

![raw](https://github.com/andreabassi78/napari-psf-simulator/raw/main/images/animation.gif)
**Slicing through a PSF aberrated with Zernike polynomials of order N=3, M=1 (coma)**

3) Click on the `Plot PSF Profile in Console` checkbox to see the x and z profiles of the PSF.
   They will show up in  the viewer console when `Calculate PSF` is executed.

![raw](https://github.com/andreabassi78/napari-psf-simulator/raw/main/images/Plot.png)
**Plot profile of the PSF, shown in the Console**

## Detailed documentation

An exhaustive documentation of the use of the plugin on scalar and vectoral propagation models can be found in [this] presentation.

A detailed explanation of the uses and advantages that simulating a PSF brings can be found [here].

The vectorial propagation model implements a secondary library: [pyfocus](https://github.com/fcaprile/PyFocus). The full documentation of this library can be found at [read the docs](https://pyfocus.readthedocs.io/en/latest/) and in the paper: "PyFocus: A Python package for vectorial calculations of focused optical fields under realistic conditions. Application to toroidal foci." https://doi.org/10.1016/j.cpc.2022.108315

## Contributing

Contributions are very welcome. Tests can be run with [tox], please ensure
the coverage at least stays the same before you submit a pull request. 
The plugin has been concived to be modular allowing the insertion of new aberations and pupils. Please contact the developers on github for adding new propagations and aberrations types. 
Any suggestions or contributions are welcome.

## License

Distributed under the terms of the [BSD-3] license,
"napari-psf-simulator" is free and open source software

## Issues

If you encounter any problems, please [file an issue] along with a detailed description.

[napari]: https://github.com/napari/napari
[Cookiecutter]: https://github.com/audreyr/cookiecutter
[@napari]: https://github.com/napari
[MIT]: http://opensource.org/licenses/MIT
[BSD-3]: http://opensource.org/licenses/BSD-3-Clause
[GNU GPL v3.0]: http://www.gnu.org/licenses/gpl-3.0.txt
[GNU LGPL v3.0]: http://www.gnu.org/licenses/lgpl-3.0.txt
[Apache Software License 2.0]: http://www.apache.org/licenses/LICENSE-2.0
[Mozilla Public License 2.0]: https://www.mozilla.org/media/MPL/2.0/index.txt
[cookiecutter-napari-plugin]: https://github.com/napari/cookiecutter-napari-plugin

[file an issue]: https://github.com/andreabassi78/napari-psf-simulator/issues

[napari]: https://github.com/napari/napari
[tox]: https://tox.readthedocs.io/en/latest/
[pip]: https://pypi.org/project/pip/
[PyPI]: https://pypi.org/

[this]: https://github.com/andreabassi78/napari-psf-simulator/raw/main/docs/napari_psf_simullator_presentation.pdf

[here]: https://github.com/andreabassi78/napari-psf-simulator/raw/main/docs/pyfocus_seminar.pptx