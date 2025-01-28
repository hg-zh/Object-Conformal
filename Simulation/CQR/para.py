import numpy as np
import random
import torch
import concurrent.futures
import os, sys
sys.path.insert(0, os.path.abspath(".."))
from cqr_comparison import ConformalizedQR, NeuralNetworkQR

import pdb
import pyreadr
import json

def compute_proportions(x_in, y_in, x_grid):
    # Initialize an output array to store the proportions
    proportions = np.zeros(len(x_grid) - 1, dtype=float)

    # Iterate through each interval defined by x_grid
    for j in range(len(x_grid) - 1):
        # Find indices of x_in that reside in the j-th interval of x_grid
        in_interval = (x_in >= x_grid[j]) & (x_in < x_grid[j + 1])

        # Compute the proportion of y_in being True for these indices
        if np.any(in_interval):  # Check if there are any x_in in this interval
            proportions[j] = np.mean(y_in[in_interval])

    return proportions

def compute_means(x_in, y_in, x_grid):
    # Initialize an output array to store the means
    means = np.zeros(len(x_grid) - 1, dtype=np.float32)

    # Iterate through each interval defined by x_grid
    for j in range(len(x_grid) - 1):
        # Find indices of x_in that reside in the j-th interval of x_grid
        in_interval = (x_in >= x_grid[j]) & (x_in < x_grid[j + 1])

        # Compute the mean of y_in for those x_in in the interval
        if np.any(in_interval):  # Only compute if there are elements in the interval
            means[j] = np.mean(y_in[in_interval])
        else:
            means[j] = np.nan  # Assign NaN if no values fall in this interval

    return means



def process_combination(setting, n, conf_methods_list, ratio_train, quantiles_net, verbose, significance, workgrid_x):
    # Load the data for the current setting and n
    Test_X = pyreadr.read_r(f"Data/Setting_heter_{setting}_n_{n}_Test_X.RData")["Testing_X"]
    Test_Y = pyreadr.read_r(f"Data/Setting_heter_{setting}_n_{n}_Test_Y.RData")["Testing_Y"]
    Train_X = pyreadr.read_r(f"Data/Setting_heter_{setting}_n_{n}_Train_X.RData")["Training_X"]
    Train_Y = pyreadr.read_r(f"Data/Setting_heter_{setting}_n_{n}_Train_Y.RData")["Training_Y"]
    
    # Convert to numpy arrays
    Train_X = np.asarray(Train_X).astype(np.float32)
    Train_Y = np.asarray(Train_Y).astype(np.float32)
    Test_X = np.asarray(Test_X).astype(np.float32)
    Test_Y = np.asarray(Test_Y).astype(np.float32)

    # Initialize coverage and length arrays
    Cov = [np.zeros((200, 100)) for _ in range(3)]
    Len = [np.zeros((200, 100)) for _ in range(3)]

    # Loop over repetitions
    for r in range(200):
        X_train = Train_X[r, :].reshape(-1, 1)
        y_train = Train_Y[r, :]
        X_test = Test_X[r, :].reshape(-1, 1)
        y_test = Test_Y[r, :]

        # Split train data
        n_train = X_train.shape[0]
        n_half = int(np.floor(n_train * ratio_train / 100.0))
        idx = np.arange(n_train)
        idx_train, idx_cal = idx[:n_half], idx[n_half:]

        # Define model parameters
        params = {
            'in_shape': X_train.shape[1],
            'epochs': 1000,
            'lr': 0.0005,
            'hidden_size': 64,
            'batch_size': 64,
            'dropout': 0.1,
            'wd': 1e-6,
            'test_ratio': 0.05,
            'random_state': r
        }

        # Train the model
        model = NeuralNetworkQR(params, quantiles_net, verbose=verbose)
        cqr = ConformalizedQR(model, model, X_train, y_train, idx_train, idx_cal, significance)

        # Compute coverage and lengths for different methods
        for i in range(3):
            conf_method = conf_methods_list[i]
            lower, upper = cqr.predict(X_test, y_test, significance, method=conf_method)
            covered = (y_test >= lower) & (y_test <= upper)
            widths = upper - lower

            Cov[i][r, :] = compute_proportions(X_test.reshape(-1), covered, workgrid_x)
            Len[i][r, :] = compute_means(X_test.reshape(-1), widths, workgrid_x)

    # Convert numpy arrays to lists for JSON serialization
    Cov_converted = [cov.tolist() for cov in Cov]
    Len_converted = [length.tolist() for length in Len]

    # Save the results to JSON
    result_dict = {"Cov": Cov_converted, "Len": Len_converted}
    output_file = f"results/Setting_heter_{setting}_n_{n}.json"
    with open(output_file, 'w') as file:
        json.dump(result_dict, file)

    return output_file

# Parallel execution setup
if __name__ == "__main__":
    settings = [2]
    ns = [1000, 2000, 4000]
    
    # Configuration (replace with your actual values)
    conf_methods_list =["CQR", "CQRm", "CQRr"] # List of your conformal methods
    ratio_train = 75  # Example value
    quantiles_net = [0.1, 0.5, 0.9]
    verbose = False
    significance = 0.1
    workgrid_x = np.linspace(-1, 1, num=101)

    # Parallel execution
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for setting in settings:
            for n in ns:
                futures.append(executor.submit(process_combination, setting, n, conf_methods_list, ratio_train, quantiles_net, verbose, significance, workgrid_x))

        for future in concurrent.futures.as_completed(futures):
            print(f"Completed: {future.result()}")




















