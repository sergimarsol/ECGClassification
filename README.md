# ECG Classification

## Overview

ECG Classification is a MATLAB-based project focused on the classification of Electrocardiogram (ECG) signals. The project extracts various features, including statistical descriptors, pNN (percentage of successive RR intervals that differ by more than 50 ms), Poincaré plots, and employs the Pan-Tompkins algorithm. The primary goals include the detection of noisy recordings and the differentiation between normal ECGs and those indicating various types of arrhythmias.

## Table of Contents

1. [Introduction](#introduction)
2. [Features](#features)
3. [Getting Started](#getting-started)
4. [MATLAB File](#matlab-file)
5. [Presentation](#presentation)
6. [Contributing](#contributing)

## Introduction

The accurate classification of ECG signals is crucial for diagnosing cardiac conditions. ECG Classification utilizes advanced signal processing techniques to analyze ECG recordings and classify them based on a variety of features. The project aims to enhance the identification of arrhythmias and improve the detection of noisy recordings.

## Features

- **Feature Extraction:** Utilizes statistical descriptors, pNN, Poincaré plots, and the Pan-Tompkins algorithm to extract relevant features from ECG signals.
- **Noise Detection:** Identifies and filters out noisy recordings to ensure accurate classification.
- **Arrhythmia Classification:** Distinguishes between normal ECGs and those indicating different types of arrhythmias using Machine Learning models and dimensionality reduction techniques such as PCA.

## Getting Started

To get started with the ECG Classification project, follow these steps:

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/your-username/ECGClassification.git
   ```

2. **MATLAB Setup:**
   Ensure you have MATLAB installed on your system. The project may require additional toolboxes or dependencies; refer to the documentation for specific requirements.

3. **Load ECG Data:**
   Prepare your ECG dataset or use the provided sample data (directories /ECG_normal and /ECG_noisy). Ensure the data is in a compatible format.

4. **Run MATLAB Script:**
   Execute the MATLAB script to perform ECG classification on the provided dataset.

## MATLAB File

The core functionality of the project is encapsulated in the MATLAB file `ECGClassification.m`. This file contains the algorithms for feature extraction, noise detection, and arrhythmia classification.

## Presentation

Explore the `ECGClassification_presentation.pdf` file to find a comprehensive presentation providing an in-depth overview of the project. The presentation covers the background, methodology and results of the project.

## Contributing

## Getting Started

If you have any questions or feedback, feel free to reach out to me.
