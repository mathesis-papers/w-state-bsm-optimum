import itertools
import json
from fractions import Fraction
from pathlib import Path
import sys

import sympy as sp

from verify import dot, norm_squared, rotation_frame


def ghz_vectors(first_frame, second_frame, sine, cosine):
    vectors = []
    for first_index, second_index in itertools.product(range(4), repeat=2):
        if first_index == second_index == 0:
            vector = (0, 0, cosine)
        elif first_index == 0:
            vector = (0, 0, second_frame[second_index - 1][2])
        elif second_index == 0:
            vector = (0, 0, first_frame[first_index - 1][2])
        else:
            left, right = first_frame[first_index - 1], second_frame[second_index - 1]
            vector = (sine * (left[0] * right[0] - left[1] * right[1]),
                      -sine * (left[0] * right[1] + left[1] * right[0]),
                      cosine * left[2] * right[2])
        vectors.append(vector)
    return vectors


def verify_ghz():
    first_amplitude, second_amplitude = sp.symbols('p r', real=True)
    ket = sp.Matrix([first_amplitude, 0, 0, 0, 0, 0, 0, second_amplitude])
    pauli = (sp.eye(2), sp.Matrix([[0, 1], [1, 0]]),
             sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1))
    norm = first_amplitude ** 2 + second_amplitude ** 2
    sine = 2 * first_amplitude * second_amplitude
    cosine = first_amplitude ** 2 - second_amplitude ** 2
    for indices in itertools.product(range(4), repeat=3):
        operator = sp.kronecker_product(*(pauli[index] for index in indices))
        actual = (ket.T.conjugate() * operator * ket)[0]
        nonidentity = [index for index in indices if index]
        expected = 0
        if not nonidentity or nonidentity == [3, 3]:
            expected = norm
        elif nonidentity in ([3], [3, 3, 3]):
            expected = cosine
        elif nonidentity == [1, 1, 1]:
            expected = sine
        elif len(nonidentity) == 3 and sorted(nonidentity) == [1, 2, 2]:
            expected = -sine
        assert sp.expand(actual - expected) == 0, indices

    left = sp.symbols('a_x a_y a_z', real=True)
    right = sp.symbols('c_x c_y c_z', real=True)
    ellipse_vector = (left[0] * right[0] - left[1] * right[1],
                      -left[0] * right[1] - left[1] * right[0],
                      left[2] * right[2])
    target = ((norm_squared(left) - left[2] ** 2) *
              (norm_squared(right) - right[2] ** 2) + left[2] ** 2 * right[2] ** 2)
    assert sp.expand(norm_squared(ellipse_vector) - target) == 0
    first_height, second_height = sp.symbols('u w', real=True)
    gap = (1 - (1 - first_height) * (1 - second_height) - first_height * second_height
           - first_height * (1 - second_height) - second_height * (1 - first_height))
    assert sp.expand(gap) == 0
    assert sp.expand(norm ** 2 - sine ** 2 - cosine ** 2) == 0

    parameters = ((Fraction(0), Fraction(1)), (Fraction(0), Fraction(-1)),
                  (Fraction(1), Fraction(0)), (Fraction(3, 5), Fraction(4, 5)),
                  (Fraction(3, 5), Fraction(-4, 5)), (Fraction(4, 5), Fraction(3, 5)),
                  (Fraction(5, 13), Fraction(12, 13)), (Fraction(12, 13), Fraction(5, 13)))
    frames = tuple(rotation_frame(quaternion) for quaternion in
                   ((1, 0, 0, 0), (1, 1, 1, 0), (2, 1, -1, 0), (3, 1, 2, -1)))
    cases = 0
    for sine_value, cosine_value in parameters:
        assert sine_value ** 2 + cosine_value ** 2 == 1
        bound = 1 + sine_value ** 2
        for first_frame, second_frame in itertools.product(frames, repeat=2):
            vectors = ghz_vectors(first_frame, second_frame, sine_value, cosine_value)
            for first_vector, second_vector in itertools.product(vectors, repeat=2):
                total = norm_squared(first_vector) + norm_squared(second_vector)
                overlap = dot(first_vector, second_vector)
                assert 2 * bound >= total
                assert bound ** 2 - bound * total + overlap ** 2 >= 0
                cases += 1
        vectors = ghz_vectors(frames[0], frames[0], sine_value, cosine_value)
        assert dot(vectors[12], vectors[5]) == 0
        assert norm_squared(vectors[12]) + norm_squared(vectors[5]) == bound
    return {'pauli_entries_checked': 64, 'polynomial_residuals': ['0', '0', '0'],
            'rational_parameter_pairs': len(parameters), 'rational_frames': len(frames),
            'ordered_parameter_frame_parity_cases': cases,
            'attainment_checked_at_each_parameter': True,
            'scope': 'Polynomial identities and the listed finite families; the all-parameter ellipsoid proof is in note.tex.'}


def verify_swapping():
    w_terms = ((0, 0, 1), (0, 1, 0), (1, 0, 0))
    outcomes = tuple(itertools.product(('Phi', 'Psi'), (1, -1)))
    bell = {('Phi', sign): {(0, 0): 1, (1, 1): sign} for sign in (1, -1)}
    bell.update({('Psi', sign): {(0, 1): 1, (1, 0): sign} for sign in (1, -1)})
    densities = {outcome: sp.zeros(4) for outcome in outcomes}
    distribution = {}
    total_probability = Fraction(0)
    for first_outcome, second_outcome in itertools.product(outcomes, repeat=2):
        entries = [0] * 4
        for first_term, second_term in itertools.product(w_terms, repeat=2):
            entries[2 * first_term[1] + second_term[1]] += (
                bell[first_outcome].get((first_term[0], second_term[0]), 0) *
                bell[second_outcome].get((first_term[2], second_term[2]), 0))
        first_kind, first_sign = first_outcome
        second_kind, second_sign = second_outcome
        if first_kind == second_kind == 'Phi':
            expected = [first_sign + second_sign, 0, 0, 1]
        elif first_kind == 'Phi':
            expected = [0, second_sign, 1, 0]
        elif second_kind == 'Phi':
            expected = [0, first_sign, 1, 0]
        else:
            expected = [first_sign + second_sign, 0, 0, 0]
        assert entries == expected
        squared_norm = sum(entry ** 2 for entry in entries)
        probability = Fraction(squared_norm, 36)
        concurrence = (Fraction(2 * abs(entries[0] * entries[3] - entries[1] * entries[2]),
                                squared_norm) if squared_norm else Fraction(0))
        distribution[concurrence] = distribution.get(concurrence, Fraction(0)) + probability
        total_probability += probability
        vector = sp.Matrix(entries)
        densities[first_outcome] += vector * vector.T / 36
    assert total_probability == 1
    assert distribution == {Fraction(4, 5): Fraction(5, 18),
                            Fraction(1): Fraction(4, 9), Fraction(0): Fraction(5, 18)}
    average = sum(concurrence * probability for concurrence, probability in distribution.items())
    assert average == Fraction(2, 3)
    spin_flip = sp.kronecker_product(sp.Matrix([[0, -sp.I], [sp.I, 0]]),
                                     sp.Matrix([[0, -sp.I], [sp.I, 0]]))
    conditional = []
    one_node_mean = sp.Integer(0)
    for (kind, sign), density in densities.items():
        expected = (sp.Matrix([[4, 0, 0, 2 * sign], [0, 2, 0, 0],
                               [0, 0, 2, 0], [2 * sign, 0, 0, 2]]) if kind == 'Phi' else
                    sp.Matrix([[4, 0, 0, 0], [0, 2, 2 * sign, 0],
                               [0, 2 * sign, 2, 0], [0, 0, 0, 0]])) / 36
        assert density == expected
        eigenvalues = (density * spin_flip * density.conjugate() * spin_flip).eigenvals()
        roots = sorted([sp.sqrt(value) for value, multiplicity in eigenvalues.items()
                        for index in range(multiplicity)], key=float, reverse=True)
        weighted_concurrence = max(sp.Integer(0), roots[0] - sum(roots[1:]))
        probability = sp.trace(density)
        concurrence = sp.simplify(weighted_concurrence / probability)
        assert concurrence == (0 if kind == 'Phi' else sp.Rational(1, 2))
        one_node_mean += weighted_concurrence
        conditional.append({'outcome': kind + ('+' if sign == 1 else '-'),
                            'probability': str(probability), 'concurrence': str(concurrence),
                            'spin_flip_square_roots_unnormalized': list(map(str, roots))})
    assert one_node_mean == sp.Rational(2, 9)
    return {'bell_outcomes_checked': 16, 'probability_sum': '1',
            'concurrence_distribution': {str(key): str(value) for key, value in sorted(distribution.items())},
            'joint_outcome_mean_concurrence': str(average),
            'one_node_outcome_mean_concurrence': str(one_node_mean),
            'one_node_outcomes': conditional}


if __name__ == '__main__':
    result = {'status': 'pass', 'arithmetic': 'exact polynomial and rational arithmetic; exact spin-flip eigenvalues',
              'generalized_ghz': verify_ghz(), 'w_swapping': verify_swapping()}
    text = json.dumps(result, indent=2, sort_keys=True) + '\n'
    print(text, end='')
    if len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(text)
