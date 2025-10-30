# Image_Processing

Code for Image Processing, finding bounding boxes and centroids.

## Overview

This repository contains MATLAB code focused on various image processing techniques, specifically targeting the identification and analysis of objects within images. The major functionalities include:

- **Finding Bounding Boxes:** Detecting objects and calculating their bounding boxes to localize them within an image.
- **Finding Centroids:** Computing the centroids of detected objects for further analysis or tracking.

## Features

- Object detection using image segmentation techniques.
- Calculation of bounding boxes for segmented objects.
- Determination of object centroids for measurement and tracking.
- Easily extensible MATLAB scripts for additional image processing tasks.

## Getting Started

### Prerequisites

- [MATLAB](https://www.mathworks.com/products/matlab.html) (R2018a or newer is recommended)

### Usage

1. Clone the repository:

    ```bash
    git clone https://github.com/syedaliabbas1/Image_Processing.git
    ```

2. Open MATLAB and add the cloned folder to your MATLAB path.

3. Run the main script (replace `main_script.m` with your entry-point file if different):

    ```matlab
    main_script
    ```

4. Modify parameters as needed for your specific image data.

### Example

```matlab
% Example usage for finding bounding boxes and centroids
img = imread('sample_image.png');
bw = imbinarize(rgb2gray(img));
props = regionprops(bw, 'BoundingBox', 'Centroid');

imshow(img); hold on;
for k = 1:length(props)
    rectangle('Position', props(k).BoundingBox, 'EdgeColor','r');
    plot(props(k).Centroid(1), props(k).Centroid(2), 'bo');
end
hold off;
```

## Repository Structure

```
Image_Processing/
├── main_script.m
├── bounding_box.m
├── centroid.m
├── sample_images/
└── README.md
```

- `main_script.m`: Main entry-point for running the image processing pipeline.
- `bounding_box.m`: Functions or scripts for bounding box detection.
- `centroid.m`: Functions or scripts for centroid calculation.
- `sample_images/`: Example images for testing.

## Contributing

Contributions are welcome! Please open issues or submit pull requests for suggestions and improvements.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Author

- [Syed Ali Abbas](https://github.com/syedaliabbas1)
