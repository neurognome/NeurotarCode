import numpy as np
from itertools import combinations


def separate_data(data, n_sections):
    section_length = round(len(data) / n_sections)
    data_sectioned = []
    for s in range(n_sections):
        if s == n_sections - 1:
            data_sectioned.append(data[s * section_length:])
        else:
            data_sectioned.append(data[s * section_length:(s + 1) * section_length])
    return data_sectioned


def bin_data(data, heading, bin_edges):
    n_bins = len(bin_edges)
    bin_counts = np.zeros(n_bins)
    binned_data = np.zeros((data.shape[1], n_bins))
    for c in range(data.shape[1]):  # Going through the cells...
        for ii in range(n_bins):
            if ii == len(bin_edges) - 1:  # Because python 0 index, so the last one is equal to length - 1
                is_currbin = heading > bin_edges[ii]  # Anything greater than the last bin, just go for it
                bin_counts[ii] = sum(is_currbin)  # Anything greater than the last bin, just go for it
            else:
                is_currbin = np.logical_and(heading > bin_edges[ii], heading < bin_edges[ii + 1])
                bin_counts[ii] = sum(is_currbin)

                binned_data[c, ii] = np.mean(data[:, c][np.squeeze(is_currbin)])
    return binned_data


class Analyzer:
    # Default class attributes
    _radius = 120  # Default cage radius

    def __init__(self, data, floating):
        self.data = dict(data)  # We call all as "dict" so that we don't affect the variables outside... make a copy
        self.floating = dict(floating)  # this is a structure...
        self._initData = dict(data)  # Initial data to recover
        self._initFloating = dict(floating)

        self._isMoving = None
        self._remove_slow_flag = False

        self.mean_quadrant_correlation = None
        self.isCorrelated = None

    def find_moving_samples(self, speed_threshold=10, peak_width=5, search_window=10):
        print('Extracting only moving frames:')
        peak_width = search_window if peak_width > search_window else peak_width
        [speed] = self.__initialize_data(self.floating, 'speed')
        self._isMoving = speed > speed_threshold
        init_moving = self._isMoving.copy()
        moving_idx = np.where(init_moving)[0]

        print('     Finding moving indices')
        for idx in moving_idx:
            pre_window = init_moving[max([0, idx - search_window]):idx]
            post_window = init_moving[(idx + 1):min([len(init_moving), idx + search_window]) + 1]  # Added the 1,
            # because we want it to INCLUDE the final member of the search window (or else it returns 9)
            self._isMoving[idx] = not ((sum(pre_window) < peak_width) and (sum(post_window) < peak_width))
        self._remove_slow_flag = True
        print('     Finished!')

    def calculate_head_direction(self, n_sections=4, bin_size=20):
        print('Calculating head direction:')
        # Initialize
        [spikes] = self.__initialize_data(self.data, 'spikes')  # The brackets here request Python to unpack the tuple
        [alpha] = self.__initialize_data(self.floating, 'alpha')
        alpha = alpha.T
        spikes_section = separate_data(spikes, n_sections)
        alpha_section = separate_data(alpha, n_sections)
        bin_edges = range(-180, 180, bin_size)

        print('     Binning neurotar data; bin size: {} degrees'.format(bin_size))
        num_cells = spikes.shape[1]
        combos = list(combinations(range(n_sections), 2))  # Get possible combinations
        binned_data = np.zeros([n_sections, len(bin_edges), num_cells])
        for s in range(n_sections):
            binned_data[s, :, :] = bin_data(spikes_section[s], alpha_section[s], bin_edges).T
        cc = np.zeros([len(combos), num_cells])  # Creating a huge list, one per cell, this is a 2D list essentially

        print('     Calculating pairwise correlations between quadrants')
        for idx, c in enumerate(combos):
            for cell in range(spikes.shape[1]):
                temp = np.corrcoef(binned_data[c[0], :, cell], binned_data[c[1], :, cell])
                cc[idx, cell] = temp[0, 1]
        self.mean_quadrant_correlation = np.mean(cc, axis=0)  # Convert to array for easy meaning across dimensions
        print('     Finished!')

    def resample_neurotar_data(self, avg_window=10):
        print('Resampling neurotar data to 2P frame rates:')
        # Initializing
        avg_window = avg_window + 1 if avg_window % 2 else avg_window  # if window is odd, make it even
        [fs, nFrames] = self.__initialize_data(self.data, 'frameRate', 'numFrames')
        if self.floating['X'].shape[0] != self._initFloating['X'].shape[0]:
            self.reset()  # in case you already ran this, it'll reset data to initial
            print('     Data reset to initial data')
        print('     Extracting time from neurotar data')
        time = self.extract_time()

        twoP_sampled_times = np.arange(1 / fs, (nFrames + 1) / fs, 1 / fs)  # Adding 1 for the total length
        neurotar_matched_indices = np.zeros([len(twoP_sampled_times)], dtype=int)
        for tt in range(len(twoP_sampled_times)):
            neurotar_matched_indices[tt] = np.argmin(
                abs(time - twoP_sampled_times[tt]))  # np.argmin returns the min idx
            # THIS wild construction is solely to match
            # python and MATLAB together. When dealing with very high precision, two mins can cause MATLAB and python to
            # disagree on the order of the mins, giving indices from MATLAB and Python that were off by 1. Didn't
            # really hurt the end result, but caused everything to be a little off. This should get around that
            # Highly recommend you don't use because it's REALLLLYY slow

            # time_diff = np.round_(abs(time - twoP_sampled_times[tt]), 5)
            # neurotar_matched_indices[tt] = min(np.where(time_diff == min(time_diff))[0])

        print('     Pulling samples from neurotar data')
        for key in self.floating:
            if key == 'time':
                continue
            curr_value = np.squeeze(self.floating[key])
            new_values = np.zeros(len(neurotar_matched_indices))
            for idx, nt in enumerate(neurotar_matched_indices):
                new_values[idx] = np.mean(  # // for floor, makes it an int
                    curr_value[np.max([0, nt - avg_window//2]):np.min([nt + avg_window//2, len(curr_value)]) + 1])
            self.floating[key] = new_values
        print('     Finished!')

    def extract_time(self):
        [time_char] = self.__initialize_data(self.floating, 'time')
        time = np.array([chr(x) for x in np.nditer(time_char)]).reshape(time_char.shape[0], time_char.shape[1])
        time = time[[3, 4, 6, 7, 9, 10, 11], :]  # this is assuming none of recording last more than an hour
        time = np.array([int(x) for x in np.nditer(time)]).reshape(time.shape[0], time.shape[1])
        # below is the ternary version of a for loop, going through the rows. transpose time so it works
        time = [np.multiply(row, [6 * 10**5, 6 * 10**4, 10**4, 10**3, 10**2, 10**1, 1]) for row in time.T]
        time = np.sum(time, axis=1) / 1000
        return time

    def reset(self):
        self.data = self._initData.copy()
        self.floating = self._initFloating.copy()

    def __initialize_data(self, datastruct, *argv):
        out = []
        for arg in argv:
            data = datastruct[arg]
            if self._remove_slow_flag:
                if len(data.shape) == 1:  # To account for weird shapes... probably a better way
                    out.append(data[self._isMoving])
                else:
                    out.append(data[self._isMoving, :])
            else:
                out.append(data)
        return out
