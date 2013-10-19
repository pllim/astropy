# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Test mags.py module."""
from __future__ import division, print_function

# THIRD-PARTY
import numpy as np

# LOCAL
from .. import mags as m
from ... import units as u
from ...tests.helper import pytest


__doctest_skip__ = ['*']

_flam = u.erg / (u.cm ** 2 * u.s * u.AA)
_fnu = u.erg / (u.cm ** 2 * u.s * u.Hz)
_flux_flam = u.Quantity([3.9135e-14, 4.0209e-14, 3.9169e-14], _flam)
_flux_fnu = u.Quantity([3.20735792e-25, 3.29903646e-25, 3.21727226e-25], _fnu)
_flux_count = u.Quantity([1214.88479883, 1248.91795446, 1217.28946691], u.count)
_flux_stmag = u.Quantity([12.41858665, 12.38919182, 12.41764379], u.mag)
_flux_abmag = u.Quantity([12.63463143, 12.60403221, 12.63128047], u.mag)
_flux_obmag = u.Quantity([-7.71133775, -7.74133477, -7.71348466], u.mag)


@pytest.mark.parametrize(
    ('value', 'ans'),
    [(20, 20),
     (_flux_obmag, _flux_obmag.value),
     (_flux_count, _flux_obmag.value)])
def test_init(value, ans):
    """Test Magnitude class initialization."""
    magobj = m.Magnitude(value)
    np.testing.assert_allclose(magobj.value, ans)
    assert magobj.unit == u.mag


@pytest.mark.parametrize(
    ('value', 'physical_unit', 'ans'),
    [(_flux_count, u.count, _flux_count),
     (_flux_stmag, u.Unit(m.ST), _flux_flam)])
def test_flux_generic(value, physical_unit, ans):
    """Test Magnitude flux method."""
    magobj = m.Magnitude(value, physical_unit)
    np.testing.assert_allclose(magobj.flux, ans.to(physical_unit))


@pytest.mark.parametrize(('value'), [_flux_flam, _flux_stmag])
def test_flux_stmag(value):
    """Test STMAG to_flux() method."""
    magobj = m.Magnitude(value, m.ST)
    np.testing.assert_allclose(magobj.flux, _flux_flam.to(m.ST))


@pytest.mark.parametrize(('value'), [_flux_fnu, _flux_abmag])
def test_toflux_abmag(value):
    """Test STMAG to_flux() method."""
    magobj = m.Magnitude(value, m.AB)
    np.testing.assert_allclose(magobj.flux, _flux_fnu.to(m.AB))
