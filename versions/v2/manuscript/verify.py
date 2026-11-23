import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp


def dot(left, right):
    return sum(first * second for first, second in zip(left, right))


def norm_squared(vector):
    return dot(vector, vector)


def triple(left, right):
    return (
        Fraction(2, 3) * (left[0] * right[2] + left[2] * right[0]),
        Fraction(2, 3) * (left[1] * right[2] + left[2] * right[1]),
        Fraction(2, 3) * (left[0] * right[0] + left[1] * right[1])
        - left[2] * right[2],
    )


def scaled_diagonal(vector):
    return tuple(Fraction(weight, 3) * entry for weight, entry in zip((2, 2, -1), vector))


def rotation_frame(quaternion):
    scalar, first, second, third = quaternion
    denominator = sum(entry * entry for entry in quaternion)
    rows = (
        (scalar * scalar + first * first - second * second - third * third,
         2 * (first * second - scalar * third), 2 * (first * third + scalar * second)),
        (2 * (first * second + scalar * third),
         scalar * scalar - first * first + second * second - third * third,
         2 * (second * third - scalar * first)),
        (2 * (first * third - scalar * second), 2 * (second * third + scalar * first),
         scalar * scalar - first * first - second * second + third * third),
    )
    return tuple(tuple(Fraction(rows[row][column], denominator) for row in range(3))
                 for column in range(3))


def correlation_vectors(first_frame, second_frame):
    vectors = []
    for first_index, second_index in itertools.product(range(4), repeat=2):
        if first_index == second_index == 0:
            vector = (Fraction(0), Fraction(0), Fraction(1, 3))
        elif first_index == 0:
            vector = scaled_diagonal(second_frame[second_index - 1])
        elif second_index == 0:
            vector = scaled_diagonal(first_frame[first_index - 1])
        else:
            vector = triple(first_frame[first_index - 1], second_frame[second_index - 1])
        vectors.append(vector)
    return vectors


def verify():
    pauli = (
        sp.eye(2), sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.diag(1, -1),
    )
    unnormalized_w = sp.Matrix([0, 1, 1, 0, 1, 0, 0, 0])
    nonzero_entries = {}
    for first_index, middle_index, last_index in itertools.product(range(4), repeat=3):
        operator = sp.kronecker_product(pauli[first_index], pauli[middle_index], pauli[last_index])
        actual = sp.simplify((unnormalized_w.T.conjugate() * operator * unnormalized_w)[0] / 3)
        indices = (first_index, middle_index, last_index)
        nonidentity = [entry for entry in indices if entry]
        if not nonidentity:
            expected = sp.Integer(1)
        elif len(nonidentity) == 1:
            expected = sp.Rational(1, 3) if nonidentity == [3] else sp.Integer(0)
        elif len(nonidentity) == 2:
            expected = (sp.Rational(2, 3) if nonidentity in ([1, 1], [2, 2])
                        else sp.Rational(-1, 3) if nonidentity == [3, 3] else sp.Integer(0))
        else:
            expected = (sp.Integer(-1) if nonidentity == [3, 3, 3]
                        else sp.Rational(2, 3) if sorted(nonidentity) in ([1, 1, 3], [2, 2, 3])
                        else sp.Integer(0))
        assert actual == expected, (indices, actual, expected)
        if actual:
            nonzero_entries[''.join(map(str, indices))] = str(actual)

    left = sp.symbols('a_x a_y a_z', real=True)
    right = sp.symbols('c_x c_y c_z', real=True)
    correlation = triple(left, right)
    bilinear = dot(left, (2 * right[0], 2 * right[1], -right[2]))
    residual = (9 * norm_squared(correlation) - bilinear ** 2
                - 4 * left[2] ** 2 - 4 * right[2] ** 2
                - 4 * right[2] ** 2 * (norm_squared(left) - 1)
                - 4 * left[2] ** 2 * (norm_squared(right) - 1))
    assert sp.expand(residual) == 0
    diagonal_residual = norm_squared((2 * left[0], 2 * left[1], -left[2])) - 4 * norm_squared(left) + 3 * left[2] ** 2
    assert sp.expand(diagonal_residual) == 0

    quaternions = (
        (1, 0, 0, 0), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1),
        (1, 1, 1, 0), (1, 1, 0, 1), (1, 0, 1, 1), (1, 1, 1, 1),
        (2, 1, -1, 0), (2, -1, 0, 1), (3, 1, 2, -1), (1, -2, 1, 2),
    )
    frames = tuple(rotation_frame(quaternion) for quaternion in quaternions)
    assert len(set(frames)) == len(frames)
    for frame in frames:
        for first_index, second_index in itertools.product(range(3), repeat=2):
            assert dot(frame[first_index], frame[second_index]) == int(first_index == second_index)
    bound = Fraction(13, 9)
    pairs_checked = 0
    equality_count = 0
    for first_frame, second_frame in itertools.product(frames, repeat=2):
        vectors = correlation_vectors(first_frame, second_frame)
        for first_vector, second_vector in itertools.product(vectors, repeat=2):
            total_squared = norm_squared(first_vector) + norm_squared(second_vector)
            overlap = dot(first_vector, second_vector)
            discriminant = total_squared ** 2 - 4 * overlap ** 2
            assert discriminant >= 0
            assert 2 * bound - total_squared >= 0
            margin = bound ** 2 - bound * total_squared + overlap ** 2
            assert margin >= 0
            pairs_checked += 1
            equality_count += int(margin == 0)

    standard = frames[0]
    vectors = correlation_vectors(standard, standard)
    assert vectors[15] == (0, 0, -1)
    assert vectors[7] == (Fraction(2, 3), 0, 0)
    assert dot(vectors[15], vectors[7]) == 0
    assert norm_squared(vectors[15]) + norm_squared(vectors[7]) == bound
    assert dot(vectors[4], vectors[8]) == 0
    assert norm_squared(vectors[4]) + norm_squared(vectors[8]) == Fraction(8, 9)
    return {
        'status': 'pass',
        'arithmetic': 'exact symbolic algebra and rational arithmetic',
        'pauli_tensor_entries_checked': 64,
        'nonzero_pauli_entries': nonzero_entries,
        'polynomial_residual': '0',
        'diagonal_norm_residual': '0',
        'rational_frames': len(frames),
        'ordered_frame_pairs': len(frames) ** 2,
        'ordered_parity_pairs_checked': pairs_checked,
        'equality_instances': equality_count,
        'attained_squared_value': '13/9',
        'single_node_attained_squared_value': '8/9',
        'finite_scope': 'The listed rational frames and all ordered central-index pairs; the all-frame theorem has the written proof in note.tex.',
    }


if __name__ == '__main__':
    result = verify()
    serialized = json.dumps(result, indent=2, sort_keys=True) + '\n'
    print(serialized, end='')
    if len(sys.argv) == 2:
        Path(sys.argv[1]).write_text(serialized)
