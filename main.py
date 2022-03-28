from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt
from typing import Union, Optional, Tuple, Iterable

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
    return result

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
    def A(self, position: Union[Tuple[float, float, float], Iterable]):
        if isinstance(position, tuple):
            r = Vec(position)
            return mu0/(4*np.pi)*(self % r.dir()) / (abs(r)*abs(r))
        elif hasattr(position, '__iter__'):
            return np.array([mu0/(4*np.pi)*self/abs(r) for r in position])
        else:

            raise Exception("An error has occurred")

    def B(self, position):
        if isinstance(position, tuple):
            r = Vec(position)
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
        rs = np.array(position)
        sizes = np.apply_along_axis(np.linalg.norm, 1, rs)
        dirs = rs / sizes[:, np.newaxis]
        return mu0 / (4 * np.pi * (sizes ** 3)[:, np.newaxis]) * (3 * dirs * (dirs @ self)[:, np.newaxis] - self)

def project_along(vector: Vec, along: Vec) -> Vec2d:
    pass

m = Dipole((0,0,1))
xs, ys = np.meshgrid(np.linspace(-1,1,60), np.linspace(-1,1,60))

xi = Vec((1,0,0))
eta = Vec((0,1,0))
zeta = cross(xi, eta)

spatial_xs = xs[:,:,np.newaxis]*xi
spatial_ys = ys[:,:,np.newaxis]*eta
rs : np.ndarray = spatial_xs + spatial_ys
rs_sizes = np.apply_along_axis(np.linalg.norm, 2, rs)

import matplotlib.pyplot as plt
shape = np.shape(spatial_xs)
B_vals = m.B(rs)
B_vals = np.apply_along_axis(np.linalg.norm, 2, B_vals)
B_vals = np.where(rs_sizes > 0.1, np.log(B_vals), -15)

fig, (ax,cbarax) = plt.subplots(1,2)
cmesh = ax.pcolormesh(xs, ys, B_vals, shading='auto')
ax.set_aspect("equal")
plt.colorbar(cmesh, cbarax, ax)
fig.show()
