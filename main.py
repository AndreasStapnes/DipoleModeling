import numpy as np
import matplotlib.pyplot as plt
from typing import Union, Optional, Tuple

class Vec(np.ndarray):
    def __new__(cls, values: Optional[Tuple[float, float, float]], **kwargs):
        obj = super().__new__(cls, shape=(3,), **kwargs)
        obj.__init__(values, **kwargs)
        return obj

    def __init__(self, values, **kwargs):
        print(values)
        super().__init__(**kwargs)
        if values is not None:
            for i in range(3): self.__setitem__(i, values[i])

    def __abs__(self):
        return np.linalg.norm(self)

    def magnitude(self):
        return abs(self)

    def dir(self):
        return self / self.magnitude()


def cross(vec1 : Vec, vec2: Vec):
    n1, = np.shape(vec1); n2, = np.shape(vec2)
    assert n1 == n2 == 3
    n = n1
    comb = np.array((vec1, vec2)).T
    result = np.ndarray(n)
    for i in range(n):
        det_obj = np.concatenate((comb[:i], comb[i+1:n]))
        result[i] = np.linalg.det(det_obj) * (-1)**i
    return result


if __name__ == "__main__":
    vec1 = Vec((1, 2, 3)); vec2 = Vec((1, 0, 1))
    print(cross(vec1, vec2))