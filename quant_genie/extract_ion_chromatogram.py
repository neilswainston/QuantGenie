'''
QuantGenie (c) University of Manchester 2018

QuantGenie is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import bisect
import sys

import pymzml
from scipy.signal import find_peaks_cwt

import matplotlib.pyplot as plt
import numpy as np


def main(args):
    '''main.'''
    handles = []

    target_mz = float(args[-2])
    error_mz = target_mz * float(args[-1]) * 1e-6

    for filename in args[:-2]:
        run = pymzml.run.Reader(filename)

        xic = [(spectrum.scan_time, _get_intensity(spectrum,
                                                   target_mz,
                                                   error_mz))
               for spectrum in run if spectrum.ms_level == 1]

        t, i = zip(*xic)

        handles.append(plt.plot(t, i, linewidth=1,
                                alpha=0.5, label=filename)[0])

        # Find and plot peaks:
        idxs = find_peaks_cwt(i, np.arange(50, 80))

        plt.scatter([t[idx] for idx in idxs],
                    [i[idx] for idx in idxs])

    plt.xlabel('RT')
    plt.ylabel('Intensity')
    plt.title('XIC: m/z=' + str(target_mz))
    plt.legend(handles=handles)
    plt.show()


def _get_intensity(spectrum, target_mz, error_mz):
    '''Gets intensity.'''
    # return interp1d(spectrum.mz, spectrum.i)(target_mz)

    lwr = target_mz - error_mz
    upp = target_mz + error_mz

    lwr_i = bisect.bisect_left(spectrum.mz, lwr)
    upp_i = bisect.bisect_right(spectrum.mz, upp, lo=lwr_i)
    intensities = spectrum.i[lwr_i:upp_i]
    return max(intensities) if len(intensities) > 0 else 0


if __name__ == '__main__':
    main(sys.argv[1:])
