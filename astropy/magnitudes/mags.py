# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""This module contains basic magnitude classes."""
from __future__ import division, print_function

# THIRD-PARTY
import numpy as np

# LOCAL
from ..extern import six
from .. import units as u


__all__ = ['Magnitude', 'MagUnit', 'ST', 'AB', 'inst', 'mag']
# 'STMag', 'ABMag'


class Magnitude(u.Quantity):
    """Class to handle generic magnitude system.

    If ``value`` given is a `~astropy.units.quantity.Quantity` but not
    in the unit of magnitude, the following conversion is done:

    .. math::

        mag = -2.5 \\times log_{10} flux

    Parameters
    ----------
    value : number, `~astropy.units.quantity.Quantity` object, or sequence of `~astropy.units.quantity.Quantity` objects.
        The numerical value of this quantity in the units given by
        unit.  If a `~astropy.units.quantity.Quantity` or sequence of them,
        creates a new `~astropy.units.quantity.Quantity` object, converting
        to magnitude as needed.

    physical_unit : Unit or Quantity
        Unit or flux corresponding to magnitude 0.  For instance,
        AB magnitude : 3630*u.Jy
        instrumental : u.count/u.second
        extinction   : u.dimensionless_unscaled (or just '') [default]

    dtype : `~numpy.dtype`, optional
        The ``dtype`` of the resulting Numpy array or scalar that will
        hold the value.  If not provided, is is determined
        automatically from the input value.

    equivalencies : list of equivalence pairs, optional
        A list of equivalence pairs. See :ref:`unit_equivalencies`.

    copy : bool, optional
        If `True` (default), then the value is copied.  Otherwise, a copy
        will only be made if :func:`__array__` returns a copy, if obj is a
        nested sequence, or if a copy is needed to satisfy ``dtype``.
        (The `False` option is intended mostly for internal use, to speed
        up initialization where it is known a copy has been made already.
        Use with care.)

    """

    _physical_unit = u.dimensionless_unscaled

    def __new__(cls, value, physical_unit=None, dtype=None,
                equivalencies=[], copy=True):

        # get input in a sane shape (merge sequences, etc.)
        value_was_quantity = isinstance(value, u.Quantity)
        value = u.Quantity(value)

        if physical_unit is None:
            if cls._physical_unit == u.dimensionless_unscaled:
                physical_unit = value.unit
            else:
                physical_unit = cls._physical_unit

        elif not isinstance(physical_unit, u.Unit):
            physical_unit = u.Unit(physical_unit)

        # if input was a Quantity, but not with magnitude units, assume it
        # is a flux and convert it to a magnitude using the physical unit
        # (this will rase UnitsError if inconsistent with physical_unit)
        if value_was_quantity and value.unit.decompose() != u.mag:
            value = -2.5 * np.log10(value / physical_unit)

        self = super(Magnitude, cls).__new__(
            cls, value.value, unit=u.mag, dtype=dtype,
            equivalencies=equivalencies, copy=copy)

        self._physical_unit = physical_unit

        return self

    @property
    def physical_unit(self):
        return self._physical_unit

    @property
    def system(self):
        return getattr(self.physical_unit, 'system', None)

    @property
    def equivalencies(self):
        return [(u.mag, self.physical_unit,
                 lambda x: 10 ** (-0.4 * x),
                 lambda x: -2.5*np.log10(x))]

    @property
    def flux(self):
        """Convert magnitude to flux.

        .. math::

            flux = 10^{-0.4 \\; mag} * physical_unit
        """
        return self.to(self.physical_unit)

    def __array_finalize__(self, obj):
        super(Magnitude, self).__array_finalize__(obj)
        if hasattr(obj, '_physical_unit'):
            self._physical_unit = obj._physical_unit

    def __quantity_view__(self, obj, unit):
        # this will need to be overriden by subclasses to preserve class
        if hasattr(unit, 'bases') and u.mag in unit.bases:
            return obj.view(Magnitude)
        else:
            return super(Magnitude, self).__quantity_view__(obj, unit)

    def __quantity_instance__(self, val, unit, **kwargs):
        # this will need to be overriden by subclasses to preserve class
        result = super(Magnitude, self).__quantity_instance__(val, unit,
                                                              **kwargs)
        if hasattr(unit, 'bases') and u.mag in unit.bases:
            result.view(Magnitude)
            result._physical_unit = self._physical_unit

        return result

    def __add__(self, other):
        # add magnitudes; this will check units of other (can be complex)
        result = super(Magnitude, self).__add__(other)
        # Make new magnitude with corresponding physical units multiplied
        result._physical_unit = (self.physical_unit *
                                 getattr(other, 'physical_unit',
                                         u.dimensionless_unscaled))
        return result

    def __radd__(self, other):
        # this cannot happen unless other is a magnitude, in which case
        # it would have reached other.__add__
        return NotImplemented

    def __iadd__(self, other):
        # add magnitudes; this will check units of other (can be complex)
        super(Magnitude, self).__iadd__(other)
        # have corresponding physical units multiplied
        self._physical_unit = (self.physical_unit *
                               getattr(other, 'physical_unit',
                                       u.dimensionless_unscaled))
        return self

    def __sub__(self, other):
        # subtract magnitudes; this will check units of other (can be complex)
        result = super(Magnitude, self).__sub__(other)
        # Make new magnitude with corresponding physical units multiplied
        result._physical_unit = (self.physical_unit /
                                 getattr(other, 'physical_unit',
                                         u.dimensionless_unscaled))
        return result

    def __rsub__(self, other):
        # this cannot happen unless other is a magnitude, in which case
        # it would have reached other.__sub__
        return NotImplemented

    def __isub__(self, other):
        # subtract magnitudes; this will check units of other (can be complex)
        super(Magnitude, self).__isub__(other)
        # have corresponding physical units divided
        self._physical_unit = (self.physical_unit /
                               getattr(other, 'physical_unit',
                                       u.dimensionless_unscaled))
        return self

    def __mul__(self, other):
        if self._physical_unit != u.dimensionless_unscaled:
            raise ValueError("Cannot multiply magnitudes of quantities which "
                             "are not dimensionless by anything")

        return super(Magnitude, self).__mul__(other)

    __rmul__ = __mul__

    def __div__(self, other):
        if self._physical_unit != u.dimensionless_unscaled:
            raise ValueError("Cannot divide magnitudes of quantities which "
                             "are not dimensionless by anything")

        return super(Magnitude, self).__div__(other)

    __truediv__ = __floordiv__ = __div__

    def __rdiv__(self, other):
        if self._physical_unit != u.dimensionless_unscaled:
            raise ValueError("Cannot use magnitudes of quantities which "
                             "are not dimensionless to divide into anything")

        return super(Magnitude, self).__div__(other)

    __rtruediv__ = __rfloordiv__ = __rdiv__

    def __eq__(self, other):
        if getattr(other, 'physical_unit', None) != self.physical_unit:
            return False

        return super(Magnitude, self).__eq__(other)

    def __ne__(self, other):
        if getattr(other, 'physical_unit', None) != self.physical_unit:
            return True

        return super(Magnitude, self).__ne__(other)

    def __gt__(self, other):
        if getattr(other, 'physical_unit', None) != self.physical_unit:
            raise u.UnitsError("Cannot compare magnitudes based on different "
                               "physical units")
        return super(Magnitude, self).__gt__(other)

    def __ge__(self, other):
        if getattr(other, 'physical_unit', None) != self.physical_unit:
            raise u.UnitsError("Cannot compare magnitudes based on different "
                               "physical units")
        return super(Magnitude, self).__ge__(other)

    def __lt__(self, other):
        if getattr(other, 'physical_unit', None) != self.physical_unit:
            raise u.UnitsError("Cannot compare magnitudes based on different "
                               "physical units")
        return super(Magnitude, self).__lt__(other)

    def __le__(self, other):
        if getattr(other, 'physical_unit', None) != self.physical_unit:
            raise u.UnitsError("Cannot compare magnitudes based on different "
                               "physical units")
        return super(Magnitude, self).__le__(other)

    def _not_implemented(self, *args, **kwargs):
        raise NotImplementedError

    dot = nansum = sum = cumsum = prod = cumprod = _not_implemented

    def __str__(self):
        return "{0} {1:s}".format(self.value,
                                  getattr(self._physical_unit, 'system', None)
                                  or self.unit.to_string())

    def __repr__(self):
        prefixstr = '<' + self.__class__.__name__ + ' '
        arrstr = np.array2string(self.view(np.ndarray), separator=',',
                                 prefix=prefixstr)
        unitstr = (getattr(self._physical_unit, 'system', None) or
                   self.unit.to_string())

        return prefixstr + arrstr + ' ' + unitstr + '>'


class MagUnit(u.CompositeUnit):

    _system = None

    def __init__(self, *args, **kwargs):
        """Define a unit corresponding to a Magnitude class"""
        system = kwargs.pop('system', None)
        # couldn't get it to work subclassing from unit, so roundabout
        unit = u.Unit(*args, **kwargs)
        super(MagUnit, self).__init__(unit.scale, unit.bases, unit.powers)
        self._system = system

    @property
    def system(self):
        return self._system

    def __mul__(self, other):
        # keep possible system attributes intact
        if isinstance(other, u.UnitBase) and other == u.dimensionless_unscaled:
            return self
        elif (self == u.dimensionless_unscaled and
              isinstance(other, MagUnit)):
            return other

        if not isinstance(other, (six.string_types, u.UnitBase)):
            try:
                return Magnitude(other, self)
            except:
                pass

        return super(MagUnit, self).__mul__(other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __div__(self, other):
        # keep possible system attributes intact
        if other == u.dimensionless_unscaled:
            return self

        return super(MagUnit, self).__div__(other)


ST = MagUnit(10.**(-0.4*21.1) * u.erg / u.cm**2 / u.s / u.AA,
             system='STmag')

AB = MagUnit(10.**(-0.4*48.60) * u.erg / u.cm**2 / u.s / u.Hz,
             system='ABmag')

inst = MagUnit(u.count / u.second, system='instmag')

mag = MagUnit(u.dimensionless_unscaled)

# class SystemMagnitude(Magnitude):
#     def __str__(self):
#         return u.Quantity.__str__(self)

#     def __repr__(self):
#         return u.Quantity.__repr__(self)


# class STMag(SystemMagnitude):
#     _physical_unit = ST


# class ABMag(SystemMagnitude):
#     _physical_unit = AB
