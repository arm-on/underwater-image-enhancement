# Underwater Image Enhancement
Table of Contents
- [Intro](#introduction)
- [Applications](#applications)
- [Method](#method)
  - [White Balancing](#white-balancing)
  - [Gamma Correction](#gamma-correction)
  - [Sharpening](#sharpening)
  - [Fusion](#fusion)
- [Usage](#usage)
- [Example](#example)
- [Persian Report](#persian-report)
- [References](#references)

# Introduction
The quality of the images taken under water may degrade due to factors such as medium scattering and absorption. The absorption substantially reduces the light energy, while the scattering causes changes in the light propagation direction. This degradation may cause the image to have low contrast or to contain a layer of fog. This repository provides a MATLAB implementation of a method to reduce the aforementioned negative effects and produce better versions of the underwater images. The method was introduced by Ancuti et al. [1].

# Applications

One may ask "Why would we even bother to enhance underwater images?". The answer is that those images can help us identify the objects existing under the water. For instance, we can detect if a cable exists somewhere and examine other underwater infrastructure. Also, knowledge of the underwater objects may be useful for marine biologists and archaeologists.

[Back to Top](#)

# Method

![Method](process.png)

## White-Balancing

To compensate for the changes in the light propagation direction caused by medium scattering, which is visible as some kind of shades of unwanted colors on an image taken under water, we apply white-balancing to the image. We chose the Gray-World Algorithm to reach this goal. This algorithm was designed based on the assumption that the color in each sensor channel averages to gray over the entire image. The formula for applying the algorithm to a given image in our case simplifies to:

$$I_{rc}(x)=I_r(x)+\alpha(\bar{I}_g-\bar{I}_r)(1-I_r(x))I_g(x)$$

where $x$ is the position of a pixel, $I_{rc}$ is the new (corrected) intensity level of that pixel corresponding to the red channel, $I_r$ and $I_g$ correspond to the red and green color channels of the image $I$. Note that the image should be normalized so that the intensity levels in each channel belong to the $[0,1]$ interval. Also, $\bar{I}_g$ and $\bar{I}_r$ represent the average intensity level of the green and red channels, respectively. The above formula was designed based on the following obversations:
1. The situation under the water does not have much effect on the green channel.
2. Based on the "Opponent Color Theory" which states that green and red are opponent colors, the green channel can be used to compensate for the distortions made to the red channel.
3. According to the Gray-World Theory, the amount of the green channel's contribution to the red channel's compensation should be proportional to the difference of red and green average intensity levels.
4. The compensation process for the red channel should only be applied to the highly distorted regions.

[Back to Top](#)

## Gamma Correction

To increase the contrast between the lighter and darker regions, we use the "Gamma Correction" procedure. This results in losing some of the details of the picture. Assuming that a pixel has an initial intensity of $i$, we will replace this intensity with $\alpha * i^{\gamma}$ where $\alpha$ and $\gamma$ are two constants.

[Back to Top](#)

## Sharpening

We prepare a sharpened version of the picture to compensate for the details lost due to the "Gamma Correction". To do this, we use the "Unsharp Masking" method. The process is done as follows: 1. Make a blurry version of the image by applying a Gaussian filter to it. 2. Calculate the difference of the original image and the blurry version. 3. Add the difference to the original image.

[Back to Top](#)

## Fusion

At this step, we combine the Gamma-Corrected and Sharpened versions of the image to obtain the enhanced image. The details of this combination is written in the paper [1].

[Back to Top](#)

# Usage

There are two folders in this repo: "code-all-in-one" and "code-section-by-section". The former contains a MATLAB script named "main_all_in_one.m" which has everything we need to give an image and produce the enhanced output. The latter one consists of many MATLAB scripts each corresponding to a part of the original code. For instance, "apply_gray_world.m" applies the Gray-World algorithm to the picture.

[Back to Top](#)

# Example

Original Image             |  Enhanced Image
:-------------------------:|:-------------------------:
![Original Image](original-sample.png)  |  ![Enhanced Image](fused-sample.png)

# Persian Report

A Persian report of the implementation is also available [here](Persian-Report.pdf).

# References
1. Ancuti, C. O., Ancuti, C., De Vleeschouwer, C., & Bekaert, P. (2017). Color balance and fusion for underwater image enhancement. IEEE Transactions on image processing, 27(1), 379-393.

[Back to Top](#)
