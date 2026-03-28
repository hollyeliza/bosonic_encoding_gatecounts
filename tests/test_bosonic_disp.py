import numpy as np
from pauli_string_formation.bosonic_disp import bosonic_disp_operator_matrix

def test_bosonic_disp_operator_matrix_matches_expected():
    d = 4
    q = bosonic_disp_operator_matrix(d)

    expected = np.array([
        [0.0, 1/np.sqrt(2), 0.0, 0.0],
        [1/np.sqrt(2), 0.0, np.sqrt(2)/np.sqrt(2), 0.0],
        [0.0, np.sqrt(2)/np.sqrt(2), 0.0, np.sqrt(3/2)],
        [0.0, 0.0, np.sqrt(3/2), 0.0],
    ])

    assert np.allclose(q, expected)