from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import Union, Optional, Tuple, Iterable, List


mu0 = 4*np.pi*1e-7

def cross(vec1: Vec, vec2: Vec):
    n1, = np.shape(vec1); n2, = np.shape(vec2)
    assert n1 == n2 == 3
    n = n1
    comb = np.array((vec1, vec2)).T
    result = np.ndarray(n)
    for i in range(n):
        det_obj = np.concatenate((comb[:i], comb[i+1:n]))
        result[i] = np.linalg.det(det_obj) * (-1)**i
    return Vec(result)

class Vec2d(np.ndarray):
    def __new__(cls, values: Optional[Tuple[float, float]], **kwargs):
        obj = super().__new__(cls, shape=(2,), **kwargs)
        obj.__init__(values, **kwargs)
        return obj
    def __init__(self, values, **kwargs):
        super().__init__(**kwargs)
        for i in range(2):
            self.__setitem__(i, values[i])
    def magnitude(self):
        return self.abs()
    def dir(self):
        return self / self.magnitude()
    def __abs__(self):
        return np.linalg.norm(self)


class Vec(np.ndarray):
    def __new__(cls, values: Optional[Tuple[float, float, float]], **kwargs):
        obj = super().__new__(cls, shape=(3,), **kwargs)
        obj.__init__(values, **kwargs)
        return obj

    def __init__(self, values, **kwargs):
        super().__init__(**kwargs)
        if values is not None:
            for i in range(3): self.__setitem__(i, values[i])

    def __abs__(self):
        return np.linalg.norm(self)

    def __mod__(self, other: Vec):
        return cross(self, other)

    def magnitude(self):
        return abs(self)

    def dir(self):
        return self / self.magnitude()


class Dipole(Vec):
    position: Vec = Vec((0,0,0))
    def A(self, position: Union[Tuple[float, float, float], Iterable]):
        if isinstance(position, tuple):
            r = Vec(position) - self.position
            return mu0/(4*np.pi)*(self % r.dir()) / (abs(r)*abs(r))
        elif hasattr(position, '__iter__'):
            return np.array([mu0/(4*np.pi)*self/abs(r-self.position) for r in position])
        else:

            raise Exception("An error has occurred")

    def B(self, position):
        if isinstance(position, tuple):
            r = Vec(position) - self.position
            r_siz = abs(r)
            return mu0/(4*np.pi*r_siz*r_siz*r_siz) * (3*r.dir()*np.dot(r.dir(), self) - self)

        elif isinstance(position, np.ndarray) and np.shape(position)[-1] == 3:
            initial_shape = np.shape(position)
            flat_position = np.reshape(position, (-1,3))
            B_flattened = self.__B_flattened__(flat_position)
            return np.reshape(B_flattened, initial_shape)

        elif hasattr(position, '__iter__'):
            return self.__B_flattened__(position)

        else:
            raise Exception("positions were of invalid type")

    def __B_flattened__(self, position):
        assert hasattr(position, '__iter__')
        rs = np.array(position) - self.position[np.newaxis, :]
        sizes = np.apply_along_axis(np.linalg.norm, 1, rs)
        dirs = rs / sizes[:, np.newaxis]
        return mu0 / (4 * np.pi * (sizes ** 3)[:, np.newaxis]) * (3 * dirs * (dirs @ self)[:, np.newaxis] - self)

class State:
    position: Vec
    velocity: Vec

    def __init__(self, position: Vec, velocity: Vec):
        self.position = position
        self.velocity = velocity

    def __add__(self, other: State):
        return State(self.position + other.position, self.velocity + other.velocity)

    def __mul__(self, other: float):
        return State(self.position*other, self.velocity*other)

    def __sub__(self, other: State):
        return self + other*(-1)

    def __truediv__(self, other: float):
        return self * (1/other)


class ChargedParticle(State):
    charge: float
    mass: float
    reacting_dipoles: List[Dipole]

    def __init__(self, position: Vec, velocity: Vec, charge: float):
        super().__init__(position, velocity)
        self.charge = charge
        self.reacting_dipoles = []

    def _force_(self):
        force = Vec((0,0,0))
        for dipole in self.reacting_dipoles:
            force += self.charge * cross(self.velocity, dipole.B(self.position))
        return force

    def _changerate_(self):
        return State(self.velocity, self._force_()/self.mass)

    def timestep(self, force: Vec, timestep=1e-3):
        pass




def RK4_step(pos, f, time: float,h: float):
    f1 = f(time, pos)
    f2 = f(time + h/2, pos+h/2*f1)
    f3 = f(time + h/2, pos+h/2*f2)
    f4 = f(time+h, pos+h*f3)
    return pos + h/6*(f1+2*f2+2*f3+f4)



def project_comp(vector: Vec, along: Vec) -> Vec:
    along = along.dir()
    return vector - (vector@along)[:,np.newaxis]*along


def most_orthogonal(vector: Vec, suggestions: List[Vec]):
    suggestions = [suggestion.dir() for suggestion in suggestions]
    return  suggestions[np.argmin(
                np.abs([np.dot(suggestion, vector) for suggestion in suggestions])
            )]


def project_along(vector: Union[Vec,np.ndarray], along: Vec, axis_tip: Optional[Vec] = None) -> Vec2d:
    xi = axis_tip
    if axis_tip is None:
        xi = most_orthogonal(along, [Vec((1, 0, 0)), Vec((0, 1, 0))])
    xi = xi.dir()
    eta = cross(along, xi).dir()
    xi_coordinates: np.ndarray = vector @ xi
    eta_coordinates: np.ndarray = vector @ eta
    twodim_crd = np.stack([xi_coordinates, eta_coordinates], axis=-1)
    twodim_crd = np.apply_along_axis(Vec2d, -1, twodim_crd)
    return twodim_crd


def spatial_span(xi: Vec, eta: Vec,
                 units_xi_min: int, units_xi_max:int, steps_xi: int,
                 units_eta_min: int, units_eta_max: int, steps_eta: int):
    xi_steps = np.linspace(units_xi_min, units_xi_max, steps_xi)
    eta_steps = np.linspace(units_eta_min, units_eta_max, steps_eta)
    xis, etas = np.meshgrid(xi_steps, eta_steps)
    rs = (np.tensordot(xis, xi, axes=0) + np.tensordot(etas, eta, axes=0))
    return rs