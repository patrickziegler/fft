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
