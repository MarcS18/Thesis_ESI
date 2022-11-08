"""
This python script performs PCA analysis on the 22 initial descriptors and outputs the principal components and
cumulative variance for each successive component.
INPUTS:
Descriptors.csv - a .csv file containing calculated descriptors for a dataset, the descriptors must be named as in the
global variable "Descs"
OUTPUTS:
Descriptors_PCA.csv - a .csv file the same as Descriptors.csv but also containing the PCA dataset
Descriptors_scree.csv - a .csv file containing the cumulative variance for each successive component (known as a scree
plot when plotted)
"""

# section 1: import modules
from sklearn import preprocessing
from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

# section 2: define inputs and outputs
dir = os.path.dirname(__file__)  # get current directory to join to files
prefix = 'TSOA_PIP'
Descriptors = pd.read_csv(os.path.join(dir, f"{prefix}_full_trimmed.csv"))  # location of descriptor file
output_PCA = os.path.join(dir, f"Descriptors_{prefix}_PCA.csv")  # location of output file for PCA
output_scree = os.path.join(dir, f"Descriptors_{prefix}_scree.csv")  # location of output file for PCA

# section 3: define PCA method

# Names of the descriptor columns in Descriptors.csv
Descs = Descriptors.columns[1:-2]
columns = [f'principal component {i + 1}' for i in range(len(Descs))]

# # Fill missing values with the mean of the descriptor
# Descriptors.fillna(Descriptors.mean(), inplace=True)


def get_PCA(Descriptors, output_PCA, ouput_scree):
    Data = Descriptors[Descs]
    # scale data
    scaler = preprocessing.StandardScaler().fit(Data)
    Data = scaler.transform(Data)
    # set up PCA with n_comp=n_desc
    print('Computing PCAs')
    pca = PCA(n_components=len(columns))
    # get components
    principalComponents = pca.fit_transform(Data)
    # make a dataframe with the PCAs in
    principalDf = pd.DataFrame(data=principalComponents, columns=columns)
    # append new data to original data
    new_data = pd.concat([Descriptors, principalDf], axis=1)
    # save data
    print('Saving PCA')
    new_data.to_csv(output_PCA, index=False)
    print('Calculating Scree Plot Data')
    # get scree plot data (cumulative variance explained by components)
    cum_scree = np.cumsum(pca.explained_variance_ratio_) * 100
    # make into dataframe and save
    print('Saving Scree Data')
    cum_scree_df = pd.DataFrame(data=[cum_scree], columns=columns)
    cum_scree_df.to_csv(output_scree, index=False)
    # Plot the data to a graph
    cum_scree_df.columns = [(i + 1) for i in range(len(columns))]
    plt.figure()
    plt.plot(cum_scree, 'o-', color='blue')
    plt.title(f'Scree Plot {prefix}')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained')
    plt.savefig(f'Scree_Plot_{prefix}.png')
    plt.show()


# section 4: run method and get PCA data
get_PCA(Descriptors, output_PCA, output_scree)
print('Finished')
