# Image-Colorization-Project
The purpose of the project is to apply transformations to certain images, changing them from grayscale tones to color using a convolutional neural network of the UNet type.

I started from an existing neural network provided by MathWorks (https://www.mathworks.com/help/images/learning-to-see-in-the-dark.html). I modified the architecture to allow training the network with pairs of images composed of grayscale images and color images on different variants such as RGB, Lab, and HSV. Lab and HSV have the advantage of separating chromaticity from brightness, unlike RGB.

The dataset used is provided by robots.ox.ac.uk ( https://www.robots.ox.ac.uk/~vgg/data/flowers/17/ ). I divided it into 1252 training images, 51 validation images, and 64 test images.

The application consists of three scripts:

Support Functions: This script is designed for data set manipulation and contains functions necessary to prepare the input and output data sets.

Images Colorization: This script loads the input data sets, along with various versions of the output data sets prepared using "Support functions". The pre-trained model is loaded in this script. The script encompasses network training configuration options, as well as the training process.

Model testing: This script is dedicated to testing the trained models and calculating the quality indicators.

Results
1. The results obtained when working simultaneously on all channels
![image](https://github.com/BigFlv/Image-Colorization-Project/assets/64215652/ee19e5ce-fdc3-4a74-9f3b-d1d02eb1abba)

Column (1) contains the desired RGB images, column (2) images generated from training to RGB, column (3) images generated from training to Lab, and column (4) images generated from training to HSV.

2. The results obtained when working separately on each channel
![image](https://github.com/BigFlv/Image-Colorization-Project/assets/64215652/8adad00f-3132-4daf-85a1-84c763f860d0)

Column (1) contains the desired RGB images, column (2) the images generated from the Lab training for 150 epochs, column (3) images generated from the Lab training for 300 epochs, column (4) images generated from HSV training over 150 epochs, and column (5) images generated from HSV training over 300 epochs.
   
3. The results obtained when working with an additional network for correction
![image](https://github.com/BigFlv/Image-Colorization-Project/assets/64215652/6b780d84-ae4a-4c00-8aa5-b68e1bfd3e29)

Column (1) contains the desired RGB images, column (2) the images generated from training to residuals correcting RGB images, column (3) corrected Lab images, column (4) corrected Lab images on channels, and column (5) corrected HSV images.
