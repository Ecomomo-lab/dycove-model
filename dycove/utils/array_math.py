"""
Functions used for cell-area weighted averaging of vegetation quantities when multiple
fractions are present in the model.
"""

import numpy as np

def get_min_max(list_of_values):
    # list_of_values: list of time-varying quantity values
    # axis 0 is list of stored time steps, axis 1 is list of mesh cells
    return np.min(list_of_values, axis=0), np.max(list_of_values, axis=0)


def cell_averaging(fractions, values):
    # sumproduct, but dividing out the magnitude of fractions because the total in each cell may be less than 1
    weighted_array = sum_product(fractions, values)
    fraction_covered = sum_elementwise(fractions)
    cell_avg = safe_array_division(weighted_array, fraction_covered)
    return cell_avg


# TODO: need to add checks for some of these functions to make sure inputs are of the correct type
def sum_product(array1, array2):    
    return sum(a*b for a, b in zip(array1, array2))

def sum_elementwise(array):
    return np.sum(array, axis=0)

def safe_array_division(array1, array2):
    array2[array2 < 0.001] = np.nan
    quotient = array1/array2
    quotient[np.isnan(array2)] = 0.
    return quotient