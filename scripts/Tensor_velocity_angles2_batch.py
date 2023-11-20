import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from openpiv import tools, validation, filters, process
from skimage.io import imread
from scipy.interpolate import RectBivariateSpline


# The program takes in microscopy images representing time series and csv
# files representing tensor fields generated from the same microscopy fields
# using Orientation J in FIJI image J. The code calculates the velocity field
# using PIV and then calculates the angles between velocity vectors and cell
# orientation tensors for each point in the grid. The data are plotted as a
# distribution of angles in a radial histogram.


def zero_to_ninety(angles):
    "convert from 0-to-180 degrees to 0-to-90 degrees"
    new_list = []
    for angle in angles:
        if angle > 90:
            angle = 180 - angle
        else:
            angle = angle
        new_list.append(angle)
    return new_list


def calculate_angle(dx, dy, u1, v1):
    # Calculate dot product and magnitudes
    dot_product = dx * u1 + dy * v1
    magnitude1 = np.sqrt(dx ** 2 + dy ** 2)
    magnitude2 = np.sqrt(u1 ** 2 + v1 ** 2)

    # Avoid division by zero
    magnitude_product = magnitude1 * magnitude2
    magnitude_product[magnitude_product == 0] = np.inf

    # Calculate cosine of angles
    cos_angles = dot_product / magnitude_product

    # Calculate angles in radians
    angles_rad = np.arccos(np.clip(cos_angles, -1.0, 1.0))

    # Convert radians to degrees
    angles_deg = np.degrees(angles_rad)

    return angles_deg


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


def convert_1d_to_2d_numpy_arrays(x, y, dx, dy):
    """Converts 1D numpy arrays to 2D numpy arrays where the number of columns and rows is dictated by the x and y variables.

    Args:
        x: A 1D numpy array of x coordinates.
        y: A 1D numpy array of y coordinates.
        dx: A 1D numpy array of dx vectors.
        dy: A 1D numpy array of dy vectors.

    Returns:
        x_2d: A 2D numpy array of x coordinates.
        y_2d: A 2D numpy array of y coordinates.
        dx_2d: A 2D numpy array of dx vectors.
        dy_2d: A 2D numpy array of dy vectors.
    """

    # Get the unique x and y values
    unique_x_values = np.unique(x)
    unique_y_values = np.unique(y)

    # Create a 2D numpy array with the unique x and y values as the axes
    x_2d, y_2d = np.meshgrid(unique_x_values, unique_y_values)

    # Reshape the dx and dy arrays to match the shape of the x_2d and y_2d arrays
    dx_2d = dx.reshape(x_2d.shape)
    dy_2d = dy.reshape(y_2d.shape)

    return x_2d, y_2d, dx_2d, dy_2d


def expand_vector_field(u, v, new_shape=(55, 55)):
    # Create a grid of coordinates for the new shape
    x_new = np.linspace(0, u.shape[1] - 1, new_shape[1])
    y_new = np.linspace(0, u.shape[0] - 1, new_shape[0])

    # Create interpolation functions for u and v using RectBivariateSpline
    spline_u = RectBivariateSpline(np.arange(u.shape[0]), np.arange(u.shape[1]), u)
    spline_v = RectBivariateSpline(np.arange(v.shape[0]), np.arange(v.shape[1]), v)

    # Evaluate the interpolation functions at the new grid of coordinates
    u_interpolated = spline_u(y_new, x_new)
    v_interpolated = spline_v(y_new, x_new)

    return u_interpolated, v_interpolated

files = ["1_1", "1_2", "1_4", "2_1", "2_2", "2_3", "2_4", "3_1", "3_2", "3_3", "3_4", "4_1", "4_2", "4_3", "4_4"]

list_of_lists = []

for file in files:
    Directory_images = "D:\\Stig\\Site" + file + "_cropped\\"
    csv_file_path = "D:\\Stig\\Site" + file + ".csv"


    names = extract_file_names_from_directory(Directory_images)
    names = names[104:106]

    window_size = 25
    overlap = 15
    search_area_size = 25

    data = pd.read_csv(csv_file_path)

    # Get unique slice values
    unique_slices = data['Slice'].unique()

    # Create directory for saving plots if it doesn't exist
    save_directory = "C:\\Users\\stigobo\\Desktop\\Data\\plots_broad\\"
    os.makedirs(save_directory, exist_ok=True)

    # Iterate through the time series

    for frame_a, frame_b, slice_value in zip(names, names[1:], unique_slices):  # Calculate PIV

        slice_data = data[data['Slice'] == slice_value]

        # Extract relevant columns and convert them to numpy arrays
        x = slice_data['X']
        x = x.to_numpy()

        y = slice_data['Y']
        y = y.to_numpy()

        dx = slice_data['DX']
        dx = dx.to_numpy()

        dy = slice_data['DY']
        dy = dy.to_numpy()

        # Concert from 1D to 2D numpy arrays
        x, y, dx, dy = convert_1d_to_2d_numpy_arrays(x, y, dx, dy)

        # Calculate the lengths of the original vectors
        original_lengths = np.sqrt(dx ** 2 + dy ** 2)

        # Normalize the vectors to a constant length
        normalized_dx = dx / original_lengths
        normalized_dy = dy / original_lengths

        # Create figure and axis

        # Perform PIV
        x1, y1, u1, v1, mag1 = piv(Directory_images + frame_a, Directory_images + frame_b,
                                   window_size=window_size,
                                   overlap=overlap,
                                   search_area_size=search_area_size)

        # adjust the shape of the velocity field to have the same shape as the tensor field
        x1, y1 = expand_vector_field(x1, y1, new_shape=(dx.shape[0], dx.shape[1]))
        u1, v1 = expand_vector_field(u1, v1, new_shape=(dx.shape[0], dx.shape[1]))

        angles = calculate_angle(dx, dy, u1, -v1)
        angles = angles.flatten()
        angles = zero_to_ninety(angles)

        angles = np.radians(angles)
        list_of_lists.append(angles)

        print("done with " + file)
        print("length " + str(len(angles)))

merged_list = np.concatenate(list_of_lists).tolist()
print("final length " +str(len(merged_list)))

plt.figure(figsize=(6, 6))
ax = plt.subplot(111, projection='polar')

# Create bins for the histogram
num_bins = 36  # Adjust the number of bins as needed
bins = np.linspace(0, 2 * np.pi, num_bins + 1)

# Plot the histogram
_, bins, _ = ax.hist(merged_list, bins=bins, color='skyblue', edgecolor='black')

# Set the direction of zero angle to be at the top
ax.set_theta_zero_location('N')

ax.set_theta_direction(-1)
ax.set_thetamin(0)
ax.set_thetamax(90)

# Set the labels and title
plt.title('28 hours', fontsize=16)
plt.ylabel('Angle count', fontsize=16)
plt.grid(True)

# Adjust font size of the numbers on the axis
ax.tick_params(axis='y', labelsize=16)  # Set font size for theta axis numbers

# Adjust the size of the degree numbers associated with each bin
ax.set_xticklabels([f'{int(np.degrees(angle))}°' for angle in ax.get_xticks()], fontsize=16)

plt.tight_layout(pad=4)

#save_path = os.path.join(save_directory, f'plot_slice_{slice_value}.png')
plt.savefig("28 hours", dpi = 300)







