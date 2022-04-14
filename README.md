# napari-psf-simulator

[![License](https://img.shields.io/pypi/l/napari-psf-simulator.svg?color=green)](https://github.com/andreabassi78/napari-psf-simulator/raw/main/LICENSE)
[![PyPI](https://img.shields.io/pypi/v/napari-psf-simulator.svg?color=green)](https://pypi.org/project/napari-psf-simulator)
[![Python Version](https://img.shields.io/pypi/pyversions/napari-psf-simulator.svg?color=green)](https://python.org)
[![tests](https://github.com/andreabassi78/napari-psf-simulator/workflows/tests/badge.svg)](https://github.com/andreabassi78/napari-psf-simulator/actions)
[![codecov](https://codecov.io/gh/andreabassi78/napari-psf-simulator/branch/main/graph/badge.svg)](https://codecov.io/gh/andreabassi78/napari-psf-simulator)
[![napari hub](https://img.shields.io/endpoint?url=https://api.napari-hub.org/shields/napari-psf-simulator)](https://napari-hub.org/plugins/napari-psf-simulator)

A plugin for the simulation of the 3D Point Spread Function of an optical systen, particularly a microscope objective.
 
Calculates the PSF as the squared Fourier Transform of the pupil and uses the angular spectrum to propagate the PSF in 3D.  
The following aberrations are included:
- phase aberration described by a Zernike polynomials with n-m coefficients
- aberration induced by a slab, with a refractive index different from the one at the object (only scalar approximation is used, polarization not considered yet).  

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


## Contributing

Contributions are very welcome. Tests can be run with [tox], please ensure
the coverage at least stays the same before you submit a pull request.

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
