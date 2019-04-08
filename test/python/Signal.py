# Copyright (C) 2019 Patrick Ziegler
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


import numpy as np


class Signal:

    def __init__(self, fs=1e3, model=None, noise=None):
        if model is None:
            model = [(1, 20), (1, 22)]
        self.fs = fs
        self.model = model
        self.noise = noise
        self.t = None
        self.x = None

    def __len__(self):
        try:
            return len(self.t)
        except TypeError:
            return -1

    def update(self, n):
        self.t = ((1 / self.fs) * np.arange(0, n))

        if self.noise is None:
            self.x = np.zeros(self.t.shape)
        else:
            self.x = self.noise[0] * np.random.randn(*self.t.shape) + self.noise[1]

        for amp, freq in self.model:
            self.x += amp * np.cos(2 * np.pi * freq * self.t)

        return self.t, self.x

    def padding(self, n):
        self.x = np.append(self.x, np.ones(n - len(self.t)))
        self.t = ((1 / self.fs) * np.arange(0, len(self.x) - 1))
        return self.t, self.x

    def spectrum(self):
        n = len(self.t)
        f = (self.fs / n) * np.arange(-n / 2, n / 2)
        X = np.fft.fftshift(np.fft.fft(self.x)) / (2 * np.sqrt(n))
        return f, X

    def show(self, f=None, X=None):
        import matplotlib.pyplot as plt

        if f is None or X is None:
            f, X = self.spectrum()

        plt.figure()

        plt.subplot(2, 1, 1)
        plt.plot(self.t, self.x)
        plt.xlabel("Time / s")
        plt.ylabel("Signal")

        plt.subplot(2, 1, 2)
        plt.plot(f, np.abs(X))
        plt.xlabel("Frequency / Hz")
        plt.ylabel("Spectrum")

        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    s1 = Signal(fs=5e3, model=[(1, 100), (1, 250), (1, 475)], noise=(1, 0))
    s1.update(4**6)
    s1.show()
