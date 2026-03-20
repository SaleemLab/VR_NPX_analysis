# -*- coding: utf-8 -*-
"""
Created on Sun Dec 14 00:09:13 2025

@author: masah
"""

from PIL import Image
import numpy as np
import os

# --- Configuration ---
# Path to your original image
input_path = r"Z:\\ibn-vision\\DATA\\SUBJECTS\\M24018\\histology\\M24018, 2nd, Left (L-M), 31.10.24.jpg"

# Output path
output_path = r"Z:\\ibn-vision\\DATA\\SUBJECTS\\M24018\\histology\\M24018, 2nd, Left (L-M), 31.10.24 (enhanced).jpg"

# --- Processing ---

# Load the image
img = Image.open(input_path)
img_array = np.array(img)

# Enhance red and blue channels
red_channel = img_array[:, :, 0]
green_channel = img_array[:, :, 1]
blue_channel = img_array[:, :, 2]

# Boost red and blue intensity
red_channel = np.clip(red_channel * 1, 0, 255)
blue_channel = np.clip(blue_channel * 2, 0, 255)

# Merge channels back together
enhanced_img_array = np.stack([red_channel, green_channel, blue_channel], axis=2).astype(np.uint8)

# Convert back to image
enhanced_img = Image.fromarray(enhanced_img_array)

# Save as PNG (lossless)
enhanced_img.save(output_path, format='PNG')

print(f"Saved enhanced image to: {output_path}")
