import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import pandas as pd

# This code plots a normalized tensor field and a velocity field overlaid on a microscopy image

def plot_normalized_tensor_field_on_image(csv_file_path, image_path):
    # Read CSV file using pandas
    data = pd.read_csv(csv_file_path)

    # Extract relevant columns
    x = data['X']
    y = data['Y']
    dx = data['DX']
    dy = data['DY']

    # Calculate the lengths of the original vectors
    original_lengths = np.sqrt(dx**2 + dy**2)

    # Normalize the vectors to a constant length
    normalized_dx = dx / original_lengths
    normalized_dy = dy / original_lengths

    # Load the microscopy image
    img = mpimg.imread(image_path)

    # Create figure and axis
    plt.figure(figsize=(10, 8))
    ax = plt.gca()

    # Display the microscopy image
    plt.imshow(np.flipud(img), extent=[x.min(), x.max(), y.min(), y.max()], cmap=("Greys"))

    # Plot normalized vectors as lines
    plt.quiver(x, y, normalized_dx, normalized_dy, angles='xy', scale_units='xy', scale=0.1, pivot = "mid", headlength=0, headaxislength=0, color='r', width=0.004, alpha=0.7)

    # Set labels and title
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Normalized Tensor Field Overlay on Microscopy Image')

    # Show plotj
    plt.show()

# Replace 'your_file_path.csv' and 'your_image_path.tif' with actual paths
csv_file_path = "C:\\Users\\stigobo\\Desktop\\Data\\test.csv"
image_path = "C:\\Users\\stigobo\\Desktop\\Data\\test.tif"
plot_normalized_tensor_field_on_image(csv_file_path, image_path)
