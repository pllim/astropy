.. _astropy_magnitudes:

*********************************
Magnitudes (`astropy.magnitudes`)
*********************************

Introduction
============

The `~astropy.magnitudes` package provides classes for representing magnitude,
as well as tools for converting between magnitude and flux.

.. warning::

    This package is currently a work-in-progress, and thus it is possible
    there will be significant API changes in later versions of Astropy.


Getting Started
===============

Define a generic magnitude object and convert it to flux:

    >>> from astropy import magnitudes as m
    >>> mag = m.Magnitude(-5)
    >>> mag
    <Magnitude -5 mag>
    >>> mag.flux
    <Quantity 100.0 >


Using `astropy.magnitudes`
==========================

Conversion between magnitude and flux, using a reference physical unit:

.. math::

    mag = -2.5 \; log_{10} (f / f_{ref})

    flux = 10^{-0.4 \; mag } f_{ref}


Generic Magnitude System
------------------------

In the generic magnitude system represented by
:class:`~astropy.magnitudes.mags.Magnitude`, the default physical unit is 1
(unscaled dimensionless).  However, if the value passed has a unit,
that unit will be used.  Alternatively, for flexibility, one can pass
on a specific physical unit, to which the value will be converted if
it has a unit already (unless that unit is already a magnitude).  One
can convert back to flux units with
:attr:`~astropy.magnitudes.mags.Magnitude.flux`, or by using the
general :meth:`~astropy.units.quantity.Quantity.to` method, for which
appropriate equivalencies are automatically set up.  Some examples::

  >>> from astropy import units as u
  >>> from astropy import magnitudes as m
  >>> mag = m.Magnitude([-5, -2.5])
  >>> mag
  <Magnitude [-5. ,-2.5] mag>
  >>> mag.flux
  <Quantity [ 100.,  10.] >
  >>> mag = m.Magnitude(100. * u.count / u.s)
  >>> mag
  <Magnitude -5.0 mag(ct / s)>
  >>> mag.flux
  <Quantity 100.0 ct / s>
  >>> mag.to(u.count / u.hr)
  <Quantity 360000.0 ct / h>

Like for quantity, a short-cut for initialisation is provided, 
:data:`~astropy.magnitudes.mags.mag`, which returns a |Magnitude|
instead of a |Quantity|:

  >>> 5. * m.mag
  <Magnitude 5.0 mag>
  >>> 5. * u.mag
  <Quantity 5.0 mag>

Magnitude Arithmetic
--------------------

Basic arithmatic with magnitudes is supported, with the physical unit
suitably adapted for addition and subtraction.  However,
multiplication and division is only allowed for magnitudes of
dimensionless quantities.  Examples::

  >>> from astropy import units as u
  >>> from astropy import magnitudes as m
  >>> mag = m.Magnitude(100. * u.count / u.s)
  >>> mag + 2.5*u.mag
  <Magnitude -2.5 mag(ct / s)>
  >>> diff = mag - m.Magnitude(10. * u.count / u.s)
  >>> diff
  <Magnitude -2.5 mag>
  >>> diff.flux
  <Quantity 10.0 >
  >>> diff * 0.1.
  <Magnitude -0.25 mag>
  >>> mag * 2.
  ValueError: Cannot multiply magnitudes of quantities which are not dimensionless by anything


STMag and ABMag
---------------

STMag and ABMag are two magnitude systems defined such that an object
with a specific constant flux distribution at all wavelengths will
have zero color at all wavelengths:

====== =================================================================== ===============
System Constant flux distribution                                          Zeropoint (mag)
====== =================================================================== ===============
STMAG  :math:`3.63 \times 10^{-9} \; erg \; cm^{-2} \; s^{-1} \; \AA^{-1}` -21.1
ABMAG  :math:`3.63 \times 10^{-20} \; erg \; cm^{-2} \; s^{-1} \; Hz^{-1}` -48.6
====== =================================================================== ===============

For these, the units system provides pre-set physical units, 
:data:`~astropy.units.astrophys.ST` and
:data:`~astropy.units.astrophys.AB`.  Using the
:class:`~astropy.magnitudes.mags.MagUnit` class, these are
converted to units that will set magnitudes,
:data:`~astropy.magnitudes.mags.ST` and :data:`~astropy.magnitudes.mags.AB`.
Examples::

  >>> from astropy import units as u
  >>> from astropy import magnitudes as m
  >>> mag = 20. * m.ST
  >>> mag
  <Magnitude 20.0 mag(ST)>
  >>> mag = m.Magnitude(u.Quantity(3.63e-9, u.erg / u.cm ** 2 / u.s / u.AA),
  ...                   m.ST)
  >>> mag                                 # doctest: +ELLIPSES
  <Magnitude 0.0002334... mag(ST)>
  >>> mag.flux                            # doctest: +ELLIPSES
  <Quantity 0.9997850193117575 ST>
  >>> mag.to(u.erg/u.cm**2/u.s/u.AA)      # doctest: +ELLIPSES
  <Quantity 3.63000000...e-09 erg / (Angstrom cm2 s)>

  >>> mag = m.Magnitude(u.Quantity(3.63e-20, u.erg / u.cm ** 2 / u.s / u.Hz),
  ...                   m.AB)
  >>> mag                                 # doctest: +ELLIPSES
  <Magnitude 0.000233... mag(AB)>
  >>> mag.flux                            # doctest: +ELLIPSES
  <Quantity 0.9997850... AB>
  >>> mag.to(u.kJy)                       # doctest: +ELLIPSES
  <Quantity 3.6299999... kJy>

Instrumental and Custom Magnitudes
----------------------------------

Another predefined magnitude is :data:`~astropy.units.magnitudes.mags.inst`::

  >>> from astropy import units as u
  >>> from astropy import magnitudes as m
  >>> mag = -10. * m.inst
  >>> mag
  <Magnitude -10.0 mag(ct / s)>
  >>> mag.flux
  <Quantity 10000.0 ct / s>

One can also define custom units for magnitudes, as follows:

  >>> mymagunit = m.MagUnit(u.photon/u.s)
  >>> mag = -5. * mymagunit
  >>> mag
  <Magnitude -5.0 mag(ph / s)>
  >>> mag.flux
  <Quantity 100. ph / s>

This can also be used to adjust zero points

  >>> zp = 0.*m.ST - (-12.*mymagunit)
  >>> zp
  <Magnitude 12.0 mag(s ST / ph)>
  >>> zp.flux
  <Quantity 1.5848931924611107e-05 s ST / ph>
  >>> mag += zp
  >>> mag
  <Magnitude 7.0 mag(ST)>
  >>> mag.flux                                 # doctest: +ELLIPSES
  <Quantity 0.0015848... ST>


See Also
========

Carroll, B. W., & Ostlie, D. A. 1996, An Introduction to Modern Astrophysics (1st ed.; Reading, MA: Addison-Wesley)

`HST ACS photometric systems <http://www.stsci.edu/hst/acs/analysis/zeropoints/#keywords>`_


Reference/API
=============

.. automodapi:: astropy.magnitudes.mags
