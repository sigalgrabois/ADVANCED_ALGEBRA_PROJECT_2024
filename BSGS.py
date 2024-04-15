from math import sqrt
from operator import mod

import numpy as np

import FiniteField
import FiniteFieldElement


def find_in_dict(dictionary, vector):
    for key, value in dictionary.items():
        for idx, vec in value:
            if np.array_equal(vec, vector):
                return key, idx
    return None, None


def create_baby_steps(l: FiniteField, g: FiniteFieldElement, m: int):
    baby_steps_dictionary = {}

    n = l.f_x_degree
    p = l.p
    base_vector = np.power(p, np.arange(n))

    for j in range(int(m)):
        result = g ** j
        vec_result = np.array(result.a).astype(int)
        # check if vec_result is in the dictionary already

        found_key_index = find_in_dict(baby_steps_dictionary, vec_result)
        if found_key_index[0] is None:
            # apply round for visibility and casting
            key = round(mod(round(np.sum(vec_result * base_vector)), int(m)))
            if key not in baby_steps_dictionary:
                baby_steps_dictionary[key] = []
            baby_steps_dictionary[key].append((j, vec_result))
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

    giant_element = g ** (-m)
    j = 0
    while j < m:
        result = h * (giant_element ** j)
        vec_result = np.array(result.a).astype(int)
        key, baby_step_index = find_in_dict(baby_steps_dictionary, vec_result)
        if key is not None and baby_step_index is not None:
            return j * m + baby_step_index
        j += 1

    raise ValueError("The result is not found.")


# Example usage:
p = 3
fx_coeff = [2, 1, 0, 0, 1]
l = FiniteField.FiniteField(p, fx_coeff)

a = [0, 1, 1, 1]
g = FiniteFieldElement.FiniteFieldElement(l, a)
l2 = FiniteField.FiniteField(p, fx_coeff)
h = FiniteFieldElement.FiniteFieldElement(l2, [0, 0, 1, 1])
result = BSGS(l, g, h)

print("The result is:")
print(result)

print("example of not found:")
a = [0, 1, 2, 2]
g = FiniteFieldElement.FiniteFieldElement(l, a)
h = FiniteFieldElement.FiniteFieldElement(l2, [0, 2, 2, 0])
try:
    result = BSGS(l, g, h)  # should raise an exception
except ValueError as e:
    print(e)
