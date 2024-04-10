from math import sqrt
from operator import mod

import numpy as np

import FiniteField
import FiniteFieldElement


# Baby stey:
# 1. create a dictinary with size m


def find_in_dict(dictionary, vector):
    # check if the vector is in the dictionary or not - the values of the dictionary are tuples where the second element is the vector
    # if the vector is in the dictionary return the first element of the tuple - the index of the vector
    # otherwise return None
    for key, value in dictionary.items():
        if value[1] == vector:
            return value[0]
        else:
            return None


def BSGS(l: FiniteField, g: FiniteFieldElement, h: FiniteFieldElement):
    """

    :param l:
    :param g:
    :param h: the result of g^x
    :return:
    """
    m = sqrt(l.field_size)
    # check that m is an integer
    if not m.is_integer():
        # throw an error if m is not an integer
        raise ValueError("The group size must be a perfect square.")
    # create a dictionary of size m
    baby_steps_dictionary = {}

    # create the vector 1, p, p^2, ..., p^(n-1)
    n = l.f_x_degree
    p = l.p
    # create a list to be the vector of 1, p, p^2, ..., p^(n-1)
    base_vector = np.array([])
    for i in range(n):
        np.append(base_vector, p ** i)

    # set the baby step from 0 to m
    for j in range(int(m)):
        result = g ** j
        vec_result = np.array(result.a)
        x = find_in_dict(baby_steps_dictionary, vec_result)
        if x is None:
            key = mod(sum(vec_result * base_vector), int(m))
            baby_steps_dictionary[key] = (j, vec_result)
        else:
            continue
    # giant element generation
    giant_element = g ** (-m)
    # loop for the giant steps
    j = 0
    while j < m - 1:
        result = h * (giant_element ** j)
        vec_result = result.a
        baby_step_index = find_in_dict(baby_steps_dictionary, vec_result)
        if baby_step_index is not None:
            return j * m + baby_step_index
        else:
            j += 1
    # throw an exception if the loop ends without finding the result
    raise ValueError("The result is not found.")


p=2
fx_coeff = [1,1,0,0,0,0,1]
l = FiniteField.FiniteField(p, fx_coeff)

a = [0,1,0,0,0,0]

g = FiniteFieldElement.FiniteFieldElement(l,a)
l2 = FiniteField.FiniteField(p,fx_coeff)
h = FiniteFieldElement.FiniteFieldElement(l2,[1,0,1,0,1,1])
result = BSGS(l, g, h)
print(result)

