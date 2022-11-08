"""
Divide a dataset into a test and training set using Scikit-learn
"""

import os
import pandas as pd
from sklearn.model_selection import train_test_split


dir = os.path.dirname(__file__)  # get current directory to join to files
prefix = 'TSOA_PIP'
Descriptors = pd.read_csv(os.path.join(dir, f"Structure_Analysis_Descriptors_{prefix}_Strict_Complete_Cleaned_Trimmed_Desc.csv"))  # location of descriptor file
output_train = os.path.join(dir, f"Descriptors_{prefix}_train.csv")  # location of output file for PCA
output_test = os.path.join(dir, f"Descriptors_{prefix}_test.csv")
output_split = os.path.join(dir, f"Structure_Analysis_Descriptors_{prefix}_Strict_Complete_Cleaned_Trimmed_Desc_Split.csv")

Descriptors = Descriptors.drop(Descriptors.columns[[-1]], axis=1)
Descriptors.set_index('Ligand', inplace=True)
# Get X (descriptors) and y (value to predict)
X = Descriptors.iloc[:, :-1]
y = Descriptors.iloc[:, -1]

# Split dataset into test and train for both X and y with a random 20% split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

# Create list of ligands in training set
train_data = []
for lig, data in Descriptors.iterrows():
    if lig in X_train.index:
        train_data.append(['Train'])
    else:
        train_data.append(['Test'])
# Add split to original dataset
train_data = pd.DataFrame(data=train_data, columns=['Set'])
Descriptors['Set'] = train_data['Set'].values
# Generate output files
X_train.to_csv(output_train)
X_test.to_csv(output_test)
Descriptors.to_csv(output_split)

print('Complete')

