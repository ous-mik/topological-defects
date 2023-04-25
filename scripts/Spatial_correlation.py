import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from openpiv import tools, validation, filters, process
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

def calculate_spatial_correlation(u, v):
    """
    Calculates the spatial correlation function for a 2D velocity field.

    Parameters:
        u (ndarray): 2D NumPy array of u-component velocity values.
        v (ndarray): 2D NumPy array of v-component velocity values.

    Returns:
        corr (ndarray): 2D NumPy array representing the spatial correlation function.
    """

    # Calculate the magnitude and angle of the velocity vectors
    mag = np.sqrt(u**2 + v**2)
    ang = np.arctan2(v, u)

    # Calculate the x and y components of the velocity fluctuations
    delta_u = u - np.mean(u)
    delta_v = v - np.mean(v)

    # Define the grid of distances to calculate the spatial correlation function
    nx, ny = u.shape
    r = np.arange(1, np.max((nx, ny))//2)
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))

    # Initialize the spatial correlation function
    corr = np.zeros(len(r))

    # Loop over all distances r and calculate the average inner product of the velocity fluctuations
    for i, rr in enumerate(r):
        for j in range(nx):
            for k in range(ny):
                if j + rr < nx:
                    du = delta_u[j, k]
                    dv = delta_v[j, k]
                    du2 = delta_u[j + rr, k]
                    dv2 = delta_v[j + rr, k]
                    inner_product = du * du2 + dv * dv2
                    corr[i] += inner_product * np.cos(ang[j, k] - ang[j + rr, k])

        corr[i] /= nx * ny * (nx - rr)

    return corr

def extract_file_names_from_directory(directory):
    """This function takes in a directory and extract the file_names contained in it"""
    file_names = []
    for subdir, dirs, files in os.walk(directory):
        for file in files:
            file_names.append(os.path.join(file))
    return sorted(file_names)

def piv(frame_a, frame_b, window_size, overlap, search_area_size):
    frame_a = tools.imread(frame_a)
    frame_b = tools.imread(frame_b)
    u, v, sig2noise = process.extended_search_area_piv(frame_a.astype(np.int32),
                                                         frame_b.astype(np.int32),
                                                         window_size=window_size,
                                                         overlap=overlap,
                                                         dt=1,
                                                         subpixel_method='gaussian',
                                                         search_area_size=search_area_size,
                                                         sig2noise_method='peak2peak')
    x, y = process.get_coordinates(image_size=frame_a.shape, window_size=window_size, overlap=overlap)
    u, v, mask = validation.sig2noise_val(u, v, sig2noise, threshold=1.05)
    u, v, mask = validation.global_val(u, v, (-5, 5), (-5, 5))
    u, v = filters.replace_outliers(u, v, method='localmean', max_iter=10, kernel_size=2)
    m = np.hypot(u, v)
    return x, y, u, v, m

def exponential_func(x, a, b, c):
    return a * np.exp(-b * x) + c

def do_it_all(path_inn, path_out, f, ps_x, ps_y, window_size, overlap, search_area_size):

    names = extract_file_names_from_directory(path_inn)
    names = names[:]  # setting the range of frames

    print(names)

    list_L = []

    for frame_a, frame_b in zip(names, names[1:]):  # Calculate PIV
        x, y, u, v, mag = piv(path_inn + frame_a, path_inn + frame_b,
                              window_size=window_size,
                              overlap=overlap,
                              search_area_size=search_area_size)

        u = np.flipud(u)
        v = np.flipud(v)

        # Calculate the spatial correlation function
        corr = calculate_spatial_correlation(u, v)

        # Normalize the correlation function
        corr_norm = corr / np.max(corr)

        # Fit an exponential curve to the normalized correlation function
        popt, pcov = curve_fit(exponential_func, np.arange(len(corr_norm)), corr_norm)

        # Generate the fitted curve
        fit_curve = exponential_func(np.arange(len(corr_norm)), *popt)

        f = interp1d(fit_curve, np.arange(len(corr_norm)))
        corr_length = f(0.53)

        scale = (ps_x / u.shape[0]) * 3.367
        corr_length = scale * corr_length
        print("new correlation length: ", corr_length)

        list_L.append(corr_length)

        print("done with " + frame_a + " " + frame_b)

    # Save list_L to a CSV file
    corr_df = pd.DataFrame(list_L)
    corr_df.to_csv(path_out + "correlation_length.csv", index=False)
    return

def main():
    path_in = "......"
    path_out = "......"
    folders = []

    for f in folders:
        do_it_all(path_in + "/" + f + "/",
                  path_out + "/Correlation_" + f + "/",
                  f,
                  ps_x=1791,
                  ps_y=1791,
                  window_size=50,
                  overlap=30,
                  search_area_size=50)

    return

if __name__ == "__main__":
    main()
