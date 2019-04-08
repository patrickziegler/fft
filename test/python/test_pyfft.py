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


import unittest
from pyfft import *
from Signal import *


class PyfftTest(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sig = Signal(fs=5e3, model=[(1, 100), (1, 250), (1, 475)], noise=(1, 0))
        self.sig.update(4**5)

    def test_rad2(self):
        f, Xr = self.sig.spectrum()
        X = fft_dif_radix2(self.sig.x)
        X = np.fft.fftshift(X) / (2 * np.sqrt(len(self.sig)))
        self.assertTrue(np.all(np.isclose(Xr, X)))

    def test_rad4(self):
        f, Xr = self.sig.spectrum()
        X = fft_dif_radix4(self.sig.x)
        X = np.fft.fftshift(X) / (2 * np.sqrt(len(self.sig)))
        self.assertTrue(np.all(np.isclose(Xr, X)))


if __name__ == '__main__':
    unittest.main()
