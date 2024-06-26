from BSGS import BSGS
from FiniteField import FiniteField
from FiniteFieldElement import FiniteFieldElement


def main():
    p = 5
    irreducible_polynom = [2, 4, 1]
    finite_field = FiniteField(p, irreducible_polynom)
    finite_field_element = FiniteFieldElement(finite_field, [4, 5])
    print("Finite field element :")
    print(finite_field_element)
    finite_field_element.pretty_printing()
    print("*********************************")
    element_order = finite_field_element.multiplicative_order()
    print(f"order of element {element_order}")
    # Find a generator (for the multiplicative group, l-{0}) element in the field
    generator_element = finite_field.find_generator()
    print(f"generator: {generator_element.poly_print()}")
    t = BSGS(finite_field, generator_element, finite_field_element)
    print(f"exponent :{t}")


if __name__ == "__main__":
    main()
