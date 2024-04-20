from math import sqrt
from operator import mod

import numpy as np

import FiniteField
import FiniteFieldElement


def find_in_dict(dictionary, vector):
    """
    This function finds a vector in a dictionary of vectors.
    :param dictionary: the dictionary to search in
    :param vector: the vector to search for
    :return: the key and index of the vector in the dictionary

    """
    for key, value in dictionary.items():
        for idx, vec in value:
            if np.array_equal(vec, vector):
                return key, idx
    return None, None


def create_baby_steps(l: FiniteField, g: FiniteFieldElement, m: int):
    """
    This function creates the baby steps dictionary for the BSGS algorithm.
    :param l: the finite field
    :param g: the generator element
    :param m: the group size
    :return: the baby steps dictionary

    """
    baby_steps_dictionary = {}  # dictionary to hold the baby steps

    n = l.f_x_degree  # extension dimension of the finite field
    p = l.p
    base_vector = np.power(p, np.arange(n))  # base vector for the field elements representation in the field extension
    for j in range(int(m)):  # iterate over the baby steps range (0, m-1)
        result = g ** j  # compute the result of the baby step
        vec_result = np.array(result.a).astype(int)
        # check if vec_result is in the dictionary already
        found_key_index = find_in_dict(baby_steps_dictionary, vec_result)
        if found_key_index[0] is None:  # if not found, add it to the dictionary
            # apply round for visibility and casting
            key = round(mod(round(np.sum(vec_result * base_vector)), int(m)))
            if key not in baby_steps_dictionary:  # if the key is not in the dictionary, add it
                baby_steps_dictionary[key] = []
            baby_steps_dictionary[key].append((j, vec_result))  # add the baby step to the dictionary
    return baby_steps_dictionary


def BSGS(l: FiniteField, g: FiniteFieldElement, h: FiniteFieldElement):
    # sanity checks for the input values
    if g.l != h.l:
        raise ValueError("The elements must be from the same field.")

    # set m to be the square root of the field size
    m = sqrt(l.field_size)
    if not m.is_integer():
        raise ValueError("The group size must be a perfect square.")

    # create the baby steps dictionary
    baby_steps_dictionary = create_baby_steps(l, g, m)

    giant_element = g ** (-m)  # compute the giant step element (g^-m)
    j = 0  # initialize the giant step index
    while j < m:
        result = h * (giant_element ** j)  # compute the result of the giant step
        vec_result = np.array(result.a).astype(int)  # convert the result to a vector
        key, baby_step_index = find_in_dict(baby_steps_dictionary,
                                            vec_result)  # find the result in the baby steps dictionary
        if key is not None and baby_step_index is not None:  # if found, return the result
            return j * m + baby_step_index  # return the exponent of the generator element that results in the given
            # element h
        j += 1  # increment the giant step index
    raise ValueError("The result is not found - are you sure the element is a generator?")
