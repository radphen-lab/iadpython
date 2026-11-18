.. |fork| image:: https://img.shields.io/badge/GitHub-RadPhen%20fork-68CA66
   :target: https://github.com/radphen-lab/iadpython
   :alt: RadPhen fork

.. |license| image:: https://img.shields.io/github/license/radphen-lab/iadpython?color=68CA66
   :target: https://github.com/radphen-lab/iadpython/blob/master/LICENSE.txt
   :alt: MIT license

.. |pypi| image:: https://img.shields.io/pypi/v/iadpython?label=upstream%20PyPI&color=68CA66
   :target: https://pypi.org/project/iadpython/
   :alt: Upstream PyPI

.. |docs| image:: https://readthedocs.org/projects/iadpython/badge
   :target: https://iadpython.readthedocs.io
   :alt: Upstream documentation

iadpython — RadPhen fork
========================

Based on `iadpython by Scott Prahl <https://github.com/scottprahl/iadpython>`_,
with forward-transport extensions maintained in the
`RadPhen fork <https://github.com/radphen-lab/iadpython>`_.

|fork| |license| |pypi| |docs|

``iadpython`` calculates forward and inverse light transport in plane-parallel
scattering materials using van de Hulst's adding-doubling method. This fork
extends the forward model to spectral multilayers with different host
refractive indices and externally supplied scattering phase functions.

The fork incorporates upstream ``main`` through commit ``33885de`` (checked
November 16, 2026), including integrating-sphere measurement formulas,
Monte Carlo lost-light iteration, inverse-search improvements, batch file
processing, and the ``iadp``, ``iadplus``, and ``iadsum`` commands. These are
shared upstream capabilities; the additions below describe what this fork
contributes beyond that integration.

The documentation and PyPI badges above refer to the upstream project.
Use the source installation below to obtain this fork's development version.

Installation
------------

Requires Python 3.10 or later. Install the current fork directly::

    python -m pip install "git+https://github.com/radphen-lab/iadpython.git@master"

For development::

    git clone https://github.com/radphen-lab/iadpython.git
    cd iadpython
    python -m pip install -e .

Runtime dependencies are NumPy, SciPy, pandas, matplotlib, and nbformat;
pip installs them automatically. ``ctypes`` is part of Python's standard library.
The optional ``iadc`` interface requires Scott Prahl's compiled
`IAD library <https://github.com/scottprahl/iad>`_. Lost-light calculations
require the external ``mc_lost`` executable; installing the Python package
alone does not provide that executable.

``pip install iadpython`` and ``conda install -c conda-forge iadpython`` select
the upstream distribution rather than this fork's source checkout.

Forward calculation
-------------------

For one homogeneous slab, specify single-scattering albedo ``a``, optical
thickness ``b``, HG anisotropy ``g``, and refractive index ``n``:

.. code-block:: python

    import iadpython as iad

    mu_s = 10.0  # scattering coefficient [1/mm]
    mu_a = 0.1   # absorption coefficient [1/mm]
    d = 1.0      # physical thickness [mm]

    a = mu_s / (mu_a + mu_s)
    b = (mu_a + mu_s) * d
    sample = iad.Sample(a=a, b=b, g=0.9, d=d, n=1.4,
                        n_above=1.5, n_below=1.5)
    ur1, ut1, uru, utu = sample.rt()
    r_ballistic, t_ballistic = sample.unscattered_rt()

``ur1`` and ``ut1`` are total hemispherical reflection/transmission for a
normally incident collimated beam. ``uru`` and ``utu`` are the corresponding
diffuse-illumination fluxes. The ballistic result contains only unscattered
light, including interface reflections and attenuation.

Fork extensions and their applications
--------------------------------------

**Scattering beyond Henyey–Greenstein.** The forward solver accepts angular
phase functions or precomputed Legendre moments. This allows externally
calculated Mie, polydisperse-particle, correlated-particle, or measured phase
functions to enter transport without replacing their higher angular moments
with a single HG parameter. The particle-scattering models themselves remain
the responsibility of the supplying library.

* ``pf_type='HG'`` uses analytic moments from ``g`` (the default).
* ``pf_type='TABULATED'`` takes a pandas DataFrame in ``pf_data``, indexed by
  scattering-angle cosine ``mu = cos(theta)``. Use unique nodes spanning
  ``[-1, 1]`` and one column per wavelength. The solver sorts the index,
  interpolates with cubic splines, integrates Legendre moments, and normalizes
  each spectrum.
* ``pf_type='MOMENTS'`` takes an array in ``pf_data``: shape ``(n_mom,)`` for
  one angular distribution, or ``(n_mom, n_wavelength)`` for spectral data.
  Supply moments proportional to ``integral(p(mu)*P_l(mu), mu=-1..1)``;
  they are normalized to ``a_0=1``. Do not multiply them by ``2*l+1``.
  For HG, the normalized moments are simply ``g**l``.

Direct moments bypass spline interpolation and moment integration inside
iadpython. This can avoid interpolation error and repeated evaluation of
expensive external scattering models. Accuracy still depends on the supplied
moments and their truncation; no general runtime speedup is claimed.

**Consistent delta-M reduction.** Custom phase functions use their actual
order-N moment for the forward-spike fraction, consistently in redistribution,
reduced albedo, and reduced optical thickness. The fraction is bounded below
unity for custom inputs to regularize the singular forward-scattering limit.
Legendre projection and redistribution assembly use array operations and
outer products. Table projection uses 32–128 Gauss–Legendre integration nodes,
selected from the requested moment count; it is not an adaptive error-controlled
integration or an exact treatment of arbitrarily sharp peaks.

**Explicit spectral multilayers.** A ``Layer`` groups one slab's ``a``, ``b``,
``g``, ``n``, ``d``, and phase function. Different layers can mix HG, tabulated,
and moment inputs. Arrays inside a layer represent wavelengths, while the
position in ``Sample(layers=[...])`` represents spatial order. This supports
dispersive hosts and different scatterers throughout a layered film or sample.
Scalars and length-one layer inputs broadcast across the sweep; incompatible
wavelength counts raise an error. Physical ``d`` is scalar per layer; transport
consumes ``b``, which must already include extinction times physical thickness.

**Internal refractive-index steps.** Unlike the common-host-index stacking
already available upstream, the fork explicitly inserts Fresnel interfaces
between layers with different real host indices. A shared Snell-invariant
angular grid maps refracted ray channels across the stack, isolates
nonpropagating channels, and includes total internal reflection and the
refractive-index scaling of radiance. This avoids interpolating transport
matrices between unrelated layer angle grids. Interface laws are explicit,
but the complete scattering solution remains a finite-quadrature approximation.

**Independent surroundings and absorbing boundaries.** ``n_outer_above`` and
``n_outer_below`` describe external media separately from the slides
``n_above`` and ``n_below``. This supports immersion and asymmetric boundary
setups. ``b_above`` and ``b_below`` specify slide optical thicknesses. Interface
power calculations accept complex indices and use Poynting-flux ratios for
transmission. A flat sample's ``n_sample_boundary`` can be supplied separately
from its real transport index; otherwise it follows changes to ``sample.n``.
These boundary extensions do not constitute a general complex-angle scattering
solver or automatically convert an extinction coefficient into bulk absorption.

**Analytic multilayer ballistic transport.** ``unscattered_rt()`` recursively
combines all slabs, interfaces, slide attenuation, and multiple internal
reflections, including wavelength sweeps. The ballistic calculation is
independent of ``quad_pts`` and provides a useful non-scattering reference.
It is an incoherent intensity calculation, without thin-film interference
or evanescent tunneling.

**Additional averages and limiting cases.** ``rt_diffuse_cone(R, T, nu_min,
nu_max)`` averages matrix-derived R/T uniformly over a cosine interval using
spline reintegration on the classic scalar-index grid. Its weighting differs
from Lambertian flux averaging; it is an approximate orientation average,
not a rough-surface scattering model. A zero-physical-thickness sample without
slides collapses to the direct interface between its surroundings. Derived
``mu_a()``, ``mu_s()``, and ``mu_sp()`` calculations also handle arrays and
mixed finite/semi-infinite optical thicknesses.

Spectral stack with mixed phase functions
-----------------------------------------

This self-contained example mixes a tabulated HG distribution, supplied
moments, and an ordinary HG layer. Replace the illustrative HG tables/moments
with those calculated by your scattering model:

.. code-block:: python

    import numpy as np
    import pandas as pd
    import iadpython as iad

    mu = np.linspace(-1.0, 1.0, 401)
    phase = pd.DataFrame({
        wavelength: (1 - g*g) / (2 * (1 + g*g - 2*g*mu)**1.5)
        for wavelength, g in [(500, 0.2), (600, 0.6)]
    }, index=mu)
    moments = 0.3 ** np.arange(129)

    layers = [
        iad.Layer(a=[0.8, 0.9], b=0.4, n=[1.3, 1.4],
                  pf_type='TABULATED', pf_data=phase),
        iad.Layer(a=0.7, b=[0.3, 0.5], n=1.6,
                  pf_type='MOMENTS', pf_data=moments),
        iad.Layer(a=0.4, b=0.2, n=1.2, g=0.1),
    ]
    stack = iad.Sample(layers=layers, quad_pts=8)
    ur1, ut1, uru, utu = stack.rt()  # each output has two wavelengths
    r_ballistic, t_ballistic = stack.unscattered_rt()

Column order defines the wavelength order; column labels do not trigger
spectral interpolation. A one-column table or one moment vector can be held
constant throughout a layered sweep. On a classic flat sample, multi-column
phase data must match other array lengths.

For a homogeneous sample, the moments path requires at least
``2*quad_pts+1`` moments. In index-mismatched stacks, each sublayer's effective
angular order may exceed the nominal ``quad_pts``, so supply additional
moments; the solver reports the required count if the array is too short.
The 129 moments above are sufficient for this example, not a universal bound.

Compatibility and inverse calculations
--------------------------------------

The merged fork retains upstream sphere measurement formulas, hot-start and
initial-simplex inputs, grid searches, iterative lost-light corrections, and
command-line tools. Its measurement helpers use the sample's ballistic solver
so boundary handling and absorbing slides follow the same conventions as the
forward model. Selecting an experimental wavelength also slices custom phase
data and layer properties and clears redistribution caches. The optional
``verbose`` flag remains available alongside upstream inverse arguments.
These integration changes do not establish joint inversion of an arbitrary
multilayer structure or MC corrections for every custom phase function.

Migration differences from upstream:

* Flat ``a``/``b``/``g`` arrays always mean wavelengths and are evaluated by
  ``rt()``. The old arrays-as-layers calls to ``rt_matrices()`` or
  ``simple_layer_matrices()`` raise an error; use ``Sample(layers=[...])``.
* A flat ``Sample`` cannot have array-valued ``n``. For index dispersion, use
  ``Layer(n=...)``. Multiwavelength stacks use ``rt()``; ``rt_matrices()`` is
  for one wavelength and materializes the layered representation in place.
* This fork's low-level ``specular_rt(..., nu, ...)`` takes ``nu`` in the
  **external incident medium**, rather than upstream's internal-slab angle.
  ``Sample.unscattered_rt()`` handles the convention for normal application use.
* ``hg_legendre(sample)`` remains available for legacy HG callers; custom
  phase functions use ``phase_legendre(sample)``.
* For infinite ``b``, the fork retains extinction-normalized coefficients:
  ``mu_a=1-a``, ``mu_s=a``, ``mu_sp=a*(1-g)``. Upstream instead uses
  scattering-normalized values (``mu_s=1``). Absolute coefficients in that
  limit must be compared using the same normalization.

See the `upstream documentation <https://iadpython.readthedocs.io>`_ for
shared functionality. The `post-merge review
<assessment/iadpython_postmerge_review.md>`_ records the current differences;
``assessment/iadpython_merge_result.md`` records the merge's preservation
and integration checks.

Numerical scope and validation
------------------------------

The full scattered solver currently supports normal incidence. General
oblique scattered transport is not implemented, even though ballistic
calculations accept oblique incidence. Transport is scalar and azimuthally
averaged; it does not propagate polarization or coherent optical phase.

Converge angular resolution and phase-function sampling for each application.
The mismatched-index grid grows with distinct index/escape thresholds, so
additional panels increase matrix size and cost. High index contrast and
strong, nearly lossless scattering can converge slowly. For the regression
stack with indices ``[1, 1.5, 1]``, ``a=1``, ``b=0.5`` and ``g=0`` in every
layer, diffuse R+T is approximately ``0.844``, ``0.913``, and ``0.953`` at
``quad_pts=8``, ``16``, and ``32``. The remaining deficit is a known numerical
limitation, not physical absorption. Arbitrary non-air diffuse normalization
and general complex-index multilayer transport also require separate validation.

The merge fixes matched-index grazing interfaces and restores the HG
compatibility entry point. Frozen pre-merge references cover mixed phase
functions, spectral index steps, ballistic transport, and asymmetric external
media. These checks establish preservation of the selected results, rather
than universal accuracy or a benchmarked speed advantage.

A focused post-merge run on November 16, 2026 passed 151 tests and two subtests,
covering layers, mixed phases, outer media, integration compatibility,
redistribution, Fresnel boundaries, thin-layer initialization, and inversion.
External C/MC calculations are not certified by that run.

License
-------

``iadpython`` is licensed under the MIT license. The upstream implementation
and scientific attribution remain credited to Scott Prahl.
