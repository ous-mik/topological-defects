import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from openpiv import tools, validation, filters, process
from skimage.io import imread
import scipy.ndimage
from scipy.spatial import distance
import matplotlib.patches

start_frame = 50
close_enough = 50
scale = 2500

def find_distance(x1, y1, x2, y2):
    a = x1, y1
    b = x2, y2
    return distance.euclidean(a, b)

def scale2(u, v, m, factor):
    u2 = np.multiply(u, factor)
    v2 = np.multiply(v, factor)
    m2 = np.multiply(m, factor)
    return u2, v2, m2

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

def angle(x, y) -> int:
    degrees = np.rad2deg(np.arctan2(y, x))
    if degrees < 0:
        degrees += 360
    return degrees

def winding_number(u, v, i, j):
    x1, y1 = u[i - 1, j - 1], v[i - 1, j - 1]
    x2, y2 = u[i,     j - 1], v[i,     j - 1]
    x3, y3 = u[i + 1, j - 1], v[i + 1, j - 1]
    x4, y4 = u[i + 1, j],     v[i + 1, j    ]
    x5, y5 = u[i + 1, j + 1], v[i + 1, j + 1]
    x6, y6 = u[i,     j + 1], v[i,     j + 1]
    x7, y7 = u[i - 1, j + 1], v[i - 1, j + 1]
    x8, y8 = u[i - 1, j    ], v[i - 1, j    ]

    x = [x1, x2, x3, x4, x5, x6, x7, x8]
    y = [y1, y2, y3, y4, y5, y6, y7, y8]

    a = []
    for i in range(8):
        a.append(angle(x[i], y[i]))

    delta = []
    for i in range(8):
        if i == 7:
            d = a[0] - a[7]
        else:
            d = a[i+1] - a[i]
        if d < -180:
            d += 360
        elif d > +180:
            d -= 360
        delta.append(d)

    skip = 0
    d_sum = 0
    if skip > 0:
        delta = sorted(delta)
    for i in range(skip, 8 - skip):
        d_sum += delta[i]
    d_sum = int(d_sum / (360 - 90 * skip))

    mean_m = np.mean(np.hypot(x, y))
    return d_sum, delta, mean_m

def cluster(k1, k2):
    n = len(k1)
    eliminated = set(())
    for d in range(n):
        for e in range(d + 1, n):
            if e not in eliminated:
                dist = find_distance(k1[d], k2[d], k1[e], k2[e])
                if dist < close_enough:
                    eliminated.add(e)

    k1_clustered = []
    k2_clustered = []
    for e in range(n):
        if e not in eliminated:
            k1_clustered.append(k1[e])
            k2_clustered.append(k2[e])

    return k1_clustered, k2_clustered

def filter(k1, k2, k1_prev, k2_prev):
    k1_filtered = []
    k2_filtered = []
    for d in range(len(k1)):
        for e in range(len(k1_prev)):
            dist = find_distance(k1[d], k2[d], k1_prev[e], k2_prev[e])
            if (dist < close_enough):
                k1_filtered.append(k1[d])
                k2_filtered.append(k2[d])
                break
    return k1_filtered, k2_filtered

def radial_elimination(list_k2, list_k1, ps_x, ps_y):
    new_list_k2, new_list_k1 = [], []
    center_x, center_y = ps_x/2, ps_y/2
    radius = ps_x/2
    for x, y in zip(list_k2, list_k1):
        if find_distance(center_x, center_y, x, y) < center_x:
            new_list_k2.append(x)
            new_list_k1.append(y)
        else:
            pass
    return new_list_k2, new_list_k1

def vorticity(u, v):
    # Calculate the vorticity as the curl of the velocity field
    dudy, dudx = np.gradient(u)
    dvdy, dvdx = np.gradient(v)
    return dvdx - dudy

def coordinatepluss(list, number):
    return np.multiply(list, number)

def do_it_all(path_inn, path_out, f, ps_x, ps_y, window_size, overlap, search_area_size):

    names = extract_file_names_from_directory(path_inn)
    names = names[:]  # setting the range of frames

    print(names)

    list_negative, list_positive = [], []
    list_mean_magnitude = []

    k1_pos_prev = []
    k2_pos_prev = []
    k1_neg_prev = []
    k2_neg_prev = []

    for frame_a, frame_b in zip(names, names[1:]):  # Calculate PIV
        x, y, u, v, mag = piv(path_inn + frame_a, path_inn + frame_b,
                              window_size=window_size,
                              overlap=overlap,
                              search_area_size=search_area_size)
        print(u.shape)

        u, v, mag = scale2(u=u, v=v, m=mag, factor=12.63)

        im = imread(path_inn + frame_a)  # Import image

        if np.isnan(u).any():  # Some director fields may contain NaN values. If this happens: notify
            print("The frames " + frame_a + " and " + frame_b + " generated NaN values and were skipped")
            pass
        else:
            # The image frame must be flipped vertically to obtain correct alignment
            im = np.flipud(im)

            # The v component of the director filed must be flipped vertically to obtain correct quiver plot
            v = np.flipud(v)
            u = np.flipud(u)

            # Calculating winding numbers at each position of an array for each time point
            # Collect the coordinates of the defects
            k1_list_negative, k2_list_negative, k1_list_positive, k2_list_positive = [], [], [], []
            rows = len(u)
            cols = len(u[1])
            for k1 in range(1, rows - 1):
                for k2 in range(1, cols - 1):
                    q, a, m = winding_number(u, v, k1, k2)
                    if q == -1:
                        # multiply image size div by array shape
                        k1_list_negative.append(k1 * (ps_x / u.shape[1]))
                        k2_list_negative.append(k2 * (ps_y / u.shape[0]))
                    elif q == 1:
                        k1_list_positive.append(k1 * (ps_x / u.shape[1]))
                        k2_list_positive.append(k2 * (ps_y / u.shape[0]))

            # Unify clusters of nearby hits into one
            k1_list_positive, k2_list_positive = cluster(k1_list_positive, k2_list_positive)
            k1_list_negative, k2_list_negative = cluster(k1_list_negative, k2_list_negative)

            # filter hits

            if (True):
                k1_pos_filtered, k2_pos_filtered = filter(k1_list_positive, k2_list_positive, k1_pos_prev, k2_pos_prev)
                k1_neg_filtered, k2_neg_filtered = filter(k1_list_negative, k2_list_negative, k1_neg_prev, k2_neg_prev)

                # save hits
                k1_pos_prev = k1_list_positive
                k2_pos_prev = k2_list_positive
                k1_neg_prev = k1_list_negative
                k2_neg_prev = k2_list_negative

                # Copy
                k1_list_positive = k1_pos_filtered
                k2_list_positive = k2_pos_filtered
                k1_list_negative = k1_neg_filtered
                k2_list_negative = k2_neg_filtered

            k2_list_positive, k1_list_positive = radial_elimination(k2_list_positive, k1_list_positive, ps_x, ps_y)
            k2_list_negative, k1_list_negative = radial_elimination(k2_list_negative, k1_list_negative, ps_x, ps_y)

            nn = 0
            k2_list_positive = np.add(k2_list_positive, nn)
            k1_list_positive = np.add(k1_list_positive, nn)
            k2_list_negative = np.add(k2_list_negative, nn)
            k1_list_negative = np.add(k1_list_negative, nn)

            w = vorticity(u, v)
            scale_factor = (ps_x/u.shape[0])*3.367
            w = w*scale_factor
            w = scipy.ndimage.gaussian_filter(w, sigma=2, order=0)

            # Scaling the vector size (sf = scaling factor)
            sf = 0.3
            rescale_x = scipy.ndimage.zoom(x, sf)
            rescale_y = scipy.ndimage.zoom(y, sf)
            rescale_u = scipy.ndimage.zoom(u, sf)
            rescale_v = scipy.ndimage.zoom(v, sf)
            rescale_mag = scipy.ndimage.zoom(mag, sf)

            # Normalizing vector size
            r = np.power(np.add(np.power(rescale_u, 2), np.power(rescale_v, 2)), 0.5)

            # Plotting data
            fig, ax = plt.subplots()
            #ax.quiver(rescale_x, rescale_y, rescale_u, rescale_v,
                #pivot="mid",
                #scale=scale,
                #color = "k",
                #alpha=1,
                #width=0.0020)
            #ax.streamplot(x, y, u, v, 5, linewidth=0.5, color="k", maxlength=4.0, minlength=0.02)
            ax.scatter(k2_list_negative, k1_list_negative, color="k", edgecolors="k", s=50)
            p = ax.imshow(im, extent=(0, ps_x, 0, ps_y), vmin=20000000, vmax=300000000, origin="lower", cmap="Greens", alpha=0.8)

            for art in ax.get_children():
                if not isinstance(art, matplotlib.patches.FancyArrowPatch):
                    continue
                art.remove()

            plt.gca().set_aspect('equal', adjustable='box')
            plt.tight_layout()
            plt.savefig(path_out + frame_a)
            plt.close()

            print("Done with " + frame_a + " and " + frame_b)

            list_negative.append(len(k2_list_negative))
            print("negative hits: " + str(len(k2_list_negative)))
            list_positive.append(len(k2_list_positive))
            print("positive hits: " + str(len(k2_list_positive)))

            mean_magnitude = np.mean(mag)
            print(mean_magnitude)
            list_mean_magnitude.append(mean_magnitude)
    combined_hits = np.add(list_negative, list_positive)
    df = pd.DataFrame([list_negative, list_positive, combined_hits, list_mean_magnitude],
                      index=['Negative', 'Positive', 'Combined', "Magnitude"])
    df2 = df.T

    df2.to_csv(path_out + "Defect_hits" + f + ".csv", index=False)

    return


def main():
    path_in = "........"
    path_out = "......."

    folders = []


    for f in folders:
        do_it_all(path_in + "/" + f + "/",
                  path_out + "\\Vorticity_" + f + "/",
                  f,
                  ps_x=1791,
                  ps_y=1791,
                  window_size=50,
                  overlap=30,
                  search_area_size=50)
    return


if __name__ == "__main__":
    main()
