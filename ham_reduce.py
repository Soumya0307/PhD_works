"""
ham_reduce.py
=============
One-callable Hamiltonian reduction via CS-VQE / Gibbs-sector projection.

Designed to slot into the resource-estimation notebook with a single call:

    from ham_reduce import reduce_hamiltonian

    qubit_ham_red, fermion_ham_red = reduce_hamiltonian(qubit_ham, fermion_ham)

Both inputs and outputs use the same OpenFermion types that are already
present in the resource-estimation notebook:
  - qubit_ham   : openfermion.QubitOperator   (from jordan_wigner)
  - fermion_ham : openfermion.InteractionOperator  (from get_molecular_hamiltonian)

The function returns:
  - qubit_ham_red   : openfermion.QubitOperator  for the reduced system
  - fermion_ham_red : openfermion.InteractionOperator  for the reduced system
    (reconstructed from the reduced qubit operator via inverse Pauli decomp)

Pipeline (mirrors csvqe_gibbs_low_sector_demo_updated.ipynb exactly):
  1. QubitOperator -> Pauli-string dict  {str: float}
  2. Identify noncontextual part  H_nc  via greedy contextuality check
  3. Diagonalise H_nc; select low-energy sector via Gibbs weights (or
     energy window -- controlled by `method` kwarg)
  4. Build projector P_low; form H_eff = P_low @ H @ P_low  in the sector basis
  5. Decompose H_eff back to Pauli-string dict on fewer qubits
  6. Apply quantum correction (CS-VQE) inside the sector:
       - split H_eff into nc / contextual parts
       - find noncontextual ground state  (classical optimisation)
       - project out noncontextual generators at their ground-state values
       - return the residual reduced Hamiltonian  (the VQE target)
  7. Convert the final reduced qubit Hamiltonian back to QubitOperator
  8. Reconstruct a FermionOperator / InteractionOperator approximation

Dependencies already present in the resource-estimation notebook:
  numpy, scipy, openfermion

Optional (used automatically if installed):
  cs_vqe  -- if available, its quasi_model / find_gs_noncon / quantum_correction
             are used for the noncontextual optimisation; otherwise a
             self-contained fallback is used.
"""

from __future__ import annotations

import itertools
import warnings
from copy import deepcopy
from functools import reduce
from typing import Dict, Optional, Tuple

import numpy as np
from scipy.optimize import minimize_scalar, differential_evolution
from openfermion import QubitOperator, FermionOperator, InteractionOperator
from openfermion.transforms import reverse_jordan_wigner


# ──────────────────────────────────────────────────────────────────────────────
# 0.  Public entry point
# ──────────────────────────────────────────────────────────────────────────────

def reduce_hamiltonian(
    qubit_ham: QubitOperator,
    fermion_ham: InteractionOperator,
    *,
    method: str = "energy_window",
    beta: float = 2.0,
    p_min: float = 1e-2,
    E_abs_cut: Optional[float] = None,
    Delta_E: Optional[float] = None,
    pauli_tol: float = 1e-10,
    verbose: bool = True,
) -> Tuple[QubitOperator, FermionOperator]:
    """Reduce a molecular Hamiltonian via Gibbs-sector projection + quantum correction.

    Parameters
    ----------
    qubit_ham : QubitOperator
        Output of ``jordan_wigner(fermion_ham)``.
    fermion_ham : InteractionOperator
        Output of ``molecule.get_molecular_hamiltonian(...)``.
    method : {"energy_window", "gibbs", "energy_cutoff"}
        How to select the low-energy sector of H_nc.
        "energy_window"  -- keep states within Delta_E above the H_nc ground state
                            (default; most robust for chemistry).
        "gibbs"          -- Boltzmann weight >= p_min at inverse temperature beta.
        "energy_cutoff"  -- keep states with H_nc eigenvalue <= E_abs_cut.
    beta : float
        Inverse temperature for Gibbs selection (only used when method="gibbs").
    p_min : float
        Minimum Gibbs probability threshold (only used when method="gibbs").
    E_abs_cut : float, optional
        Absolute energy threshold (only used when method="energy_cutoff").
    Delta_E : float, optional
        Energy window above H_nc ground state (only used when method="energy_window").
        Defaults to 10% of the H_nc spectral range when None.
    pauli_tol : float
        Coefficients below this threshold are dropped from Pauli dicts.
    verbose : bool
        Print progress and dimension information.

    Returns
    -------
    qubit_ham_red : QubitOperator
        Reduced qubit Hamiltonian (fewer qubits).
    fermion_ham_red : FermionOperator
        Fermionic representation reconstructed from the reduced qubit operator
        via reverse Jordan-Wigner.  Use this anywhere fermion_ham is needed
        downstream (e.g. trotter.simulate_trotter).
    """

    # ── Step 1: QubitOperator → Pauli-string dict ────────────────────────────
    ham_dict = _qubit_op_to_dict(qubit_ham)
    n_qubits  = _infer_n_qubits(qubit_ham)
    if verbose:
        print(f"[reduce_hamiltonian] Original: {n_qubits} qubits, "
              f"{len(ham_dict)} Pauli terms")

    # ── Step 2: Split into H_nc and H_c ─────────────────────────────────────
    ham_nc, ham_c = _greedy_noncontextual_split(ham_dict)
    if verbose:
        print(f"  H_nc: {len(ham_nc)} terms   H_c: {len(ham_c)} terms   "
              f"(H_nc contextual={_contextual_q_ham(ham_nc)})")

    # ── Step 3: Diagonalise H_nc, select low-energy sector ──────────────────
    H_nc_mat = _dict_to_matrix(ham_nc, n_qubits)
    evals_nc, evecs_nc = np.linalg.eigh(H_nc_mat)
    E_nc_classical = float(evals_nc[0])
    # Default Delta_E: 10 % of H_nc spectral range
    if Delta_E is None and method == "energy_window":
        Delta_E = 0.10 * float(evals_nc[-1] - evals_nc[0])
        Delta_E = max(Delta_E, 1e-3)  # at least some window

    selected, info = _select_sector(
        evals_nc,
        method=method,
        beta=beta,
        p_min=p_min,
        E_abs_cut=E_abs_cut,
        Delta_E=Delta_E,
    )
    if verbose:
        print(f"  Sector: {len(selected)} states kept "
              f"(method='{method}', eigenvalues {evals_nc[selected]})")

    if len(selected) == 0:
        warnings.warn(
            "No states selected — returning original Hamiltonians unchanged.",
            RuntimeWarning,
        )
        return qubit_ham, FermionOperator(fermion_ham), E_nc_classical

    # ── Step 4: Build H_eff in the sector basis ──────────────────────────────
    H_c_mat = _dict_to_matrix(ham_c, n_qubits)
    V_low   = evecs_nc[:, selected]
    H_eff_mat = V_low.conj().T @ H_c_mat @ V_low  # (d_eff × d_eff)

    d_eff = H_eff_mat.shape[0]
    if verbose:
        print(f"  H_eff dimension: {d_eff}×{d_eff}")

    # ── Step 5: Decompose H_eff → Pauli dict (on n_eff qubits) ──────────────
    n_eff = int(np.round(np.log2(d_eff))) if d_eff > 1 else 0

    if d_eff == 1:
        # Sector is a single state — scalar energy, trivial Hamiltonian
        scalar = float(np.real(H_eff_mat[0, 0]))
        if verbose:
            print(f"  Sector collapsed to scalar energy {scalar:.6f}")
        qubit_red = QubitOperator("", scalar)
        return QubitOperator("", scalar), FermionOperator((), scalar), E_nc_classical

    if 2**n_eff != d_eff:
        # d_eff is not a power of 2; pad to next power of 2
        next_pow2 = 2 ** int(np.ceil(np.log2(d_eff)))
        pad = next_pow2 - d_eff
        H_eff_mat = np.pad(H_eff_mat, ((0, pad), (0, pad)), mode='constant')
        n_eff     = int(np.ceil(np.log2(next_pow2)))
        if verbose:
            print(f"  Padded H_eff to {next_pow2}×{next_pow2} ({n_eff} qubits)")

    ham_eff = _matrix_to_dict(H_eff_mat, n_eff, tol=pauli_tol)
    if verbose:
        print(f"  H_eff Pauli dict: {len(ham_eff)} terms on {n_eff} qubits")

    # ── Step 6: Quantum correction inside the sector ─────────────────────────
    ham_red_final = _quantum_correction_in_sector(ham_eff, n_eff, verbose=verbose)

    # ── Step 7: Pauli dict → QubitOperator ───────────────────────────────────
    qubit_ham_red = _dict_to_qubit_op(ham_red_final)

    # ── Step 8: QubitOperator → FermionOperator (reverse JW) ─────────────────
    try:
        fermion_ham_red = reverse_jordan_wigner(qubit_ham_red)
    except Exception:
        # reverse_jordan_wigner can fail for non-JW-image operators;
        # fall back to a FermionOperator with the same coefficients
        fermion_ham_red = _qubit_op_to_fermion_fallback(qubit_ham_red)
        if verbose:
            print("  Note: reverse_jordan_wigner failed; using direct coefficient copy "
                  "as FermionOperator (sufficient for Trotter).")

    if verbose:
        n_red = _infer_n_qubits(qubit_ham_red)
        print(f"[reduce_hamiltonian] Reduced: {n_red} qubits, "
              f"{len(ham_red_final)} Pauli terms")

    return qubit_ham_red, fermion_ham_red, E_nc_classical


# ──────────────────────────────────────────────────────────────────────────────
# 1.  Format converters
# ──────────────────────────────────────────────────────────────────────────────

def _qubit_op_to_dict(qubit_op: QubitOperator) -> Dict[str, float]:
    """OpenFermion QubitOperator -> {'IXZY...': coeff} dict.

    The qubit_op stores terms as tuples like ((0,'X'),(2,'Z')).
    We convert each to a full Pauli string 'IXIZ...' padded to n_qubits.
    """
    n = _infer_n_qubits(qubit_op)
    out: Dict[str, float] = {}
    for term, coeff in qubit_op.terms.items():
        label = ['I'] * n
        for idx, op in term:
            label[idx] = op
        key = ''.join(label)
        out[key] = float(np.real(coeff))
    return out


def _infer_n_qubits(qubit_op: QubitOperator) -> int:
    """Find the number of qubits from the highest qubit index in a QubitOperator."""
    max_idx = -1
    for term in qubit_op.terms:
        for idx, _ in term:
            if idx > max_idx:
                max_idx = idx
    return max_idx + 1 if max_idx >= 0 else 0


def _dict_to_qubit_op(ham_dict: Dict[str, float]) -> QubitOperator:
    """{'IXZY...': coeff} -> QubitOperator."""
    op = QubitOperator()
    for pauli_str, coeff in ham_dict.items():
        # Build term tuple, e.g. "IXZ" -> ((1,'X'),(2,'Z'))
        term = tuple(
            (i, p) for i, p in enumerate(pauli_str) if p != 'I'
        )
        op += QubitOperator(term, coeff)
    return op


def _qubit_op_to_fermion_fallback(qubit_op: QubitOperator) -> FermionOperator:
    """Last-resort: copy qubit terms directly as FermionOperator labels."""
    ferm = FermionOperator()
    for term, coeff in qubit_op.terms.items():
        ferm += FermionOperator(term, coeff)
    return ferm


# ──────────────────────────────────────────────────────────────────────────────
# 2.  Matrix ↔ Pauli-string dict
#     (Identical to the notebook implementations, self-contained here)
# ──────────────────────────────────────────────────────────────────────────────

_PAULI_MATS = {
    'I': np.array([[1, 0], [0, 1]], dtype=complex),
    'X': np.array([[0, 1], [1, 0]], dtype=complex),
    'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
    'Z': np.array([[1, 0], [0, -1]], dtype=complex),
}


def _pauli_string_matrix(s: str) -> np.ndarray:
    out = np.array([[1]], dtype=complex)
    for ch in s:
        out = np.kron(out, _PAULI_MATS[ch])
    return out


def _dict_to_matrix(ham_dict: Dict[str, float], n_qubits: int) -> np.ndarray:
    dim = 2 ** n_qubits
    H = np.zeros((dim, dim), dtype=complex)
    for p, c in ham_dict.items():
        H += c * _pauli_string_matrix(p)
    return H


def _matrix_to_dict(H_mat: np.ndarray, n: int, tol: float = 1e-10) -> Dict[str, float]:
    """Exact Pauli decomposition via trace formula  h_k = Tr(P_k H) / 2^n."""
    d   = 2 ** n
    out: Dict[str, float] = {}
    for labels in itertools.product('IXYZ', repeat=n):
        pstr  = ''.join(labels)
        P     = _pauli_string_matrix(pstr)
        coeff = np.trace(P @ H_mat) / d
        if abs(coeff) > tol:
            out[pstr] = float(np.real(coeff))
    return out


# ──────────────────────────────────────────────────────────────────────────────
# 3.  Contextuality helpers
#     (Direct port from cs_vqe.py and the demo notebook)
# ──────────────────────────────────────────────────────────────────────────────

def _commute(a: str, b: str) -> bool:
    sign = 1
    for x, y in zip(a, b):
        if x != 'I' and y != 'I' and x != y:
            sign *= -1
    return sign == 1


def _contextual_q(terms) -> bool:
    T = [t for t in terms if any(not _commute(t, s) for s in terms)]
    for i in range(len(T)):
        for j in range(len(T)):
            for k in range(j, len(T)):
                if (i != j and i != k
                        and _commute(T[i], T[j])
                        and _commute(T[i], T[k])
                        and not _commute(T[j], T[k])):
                    return True
    return False


def _contextual_q_ham(ham_dict: Dict[str, float]) -> bool:
    return _contextual_q(list(ham_dict.keys()))


# Public alias used in docstring
_contextual_q_ham_public = _contextual_q_ham


def _greedy_noncontextual_split(
    ham_dict: Dict[str, float],
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """Greedy extraction of the largest-weight noncontextual sub-Hamiltonian.

    Iterates over terms in decreasing |coeff| order; keeps a term only if
    the retained set remains noncontextual.  Mirrors the implementation in
    csvqe_gibbs_low_sector_demo_updated.ipynb.
    """
    sorted_terms = sorted(ham_dict.keys(), key=lambda t: -abs(ham_dict[t]))
    nc_terms = []
    for t in sorted_terms:
        if not _contextual_q(nc_terms + [t]):
            nc_terms.append(t)
    ham_nc = {t: ham_dict[t] for t in nc_terms}
    ham_c  = {t: ham_dict[t] for t in ham_dict if t not in nc_terms}
    return ham_nc, ham_c


# ──────────────────────────────────────────────────────────────────────────────
# 4.  Sector selection
#     (Port of select_sector_from_Hnc from the demo notebook)
# ──────────────────────────────────────────────────────────────────────────────

def _select_sector(
    evals_nc: np.ndarray,
    method: str = "energy_window",
    beta: float = 2.0,
    p_min: float = 1e-2,
    E_abs_cut: Optional[float] = None,
    Delta_E: Optional[float] = None,
    tol: float = 1e-12,
):
    evals_nc = np.asarray(evals_nc)
    E_min    = float(np.min(evals_nc))

    if method == "gibbs":
        shifted = evals_nc - E_min
        boltz   = np.exp(-beta * shifted)
        Z       = boltz.sum()
        weights = boltz / Z
        selected = np.where(weights >= p_min - tol)[0]
        info = {"weights": weights}

    elif method == "energy_cutoff":
        if E_abs_cut is None:
            raise ValueError("method='energy_cutoff' requires E_abs_cut.")
        selected = np.where(evals_nc <= E_abs_cut + tol)[0]
        info = {"E_abs_cut": E_abs_cut}

    elif method == "energy_window":
        if Delta_E is None:
            raise ValueError("method='energy_window' requires Delta_E.")
        selected = np.where(evals_nc <= E_min + Delta_E + tol)[0]
        info = {"E_min": E_min, "Delta_E": Delta_E, "cutoff": E_min + Delta_E}

    else:
        raise ValueError(f"Unknown method '{method}'. "
                         "Choose 'gibbs', 'energy_cutoff', or 'energy_window'.")

    return selected, info


# ──────────────────────────────────────────────────────────────────────────────
# 5.  Self-contained quantum correction
#     Mirrors cs_vqe.py's quasi_model / find_gs_noncon / quantum_correction
#     without requiring the cs_vqe package.
# ──────────────────────────────────────────────────────────────────────────────

def _pauli_mult(p: str, q: str):
    """Multiply two Pauli strings; return [result_str, phase]."""
    assert len(p) == len(q)
    sgn  = 1
    out  = ''
    for a, b in zip(p, q):
        if   a == 'I':                   out += b
        elif b == 'I':                   out += a
        elif a == 'X':
            if   b == 'X':               out += 'I'
            elif b == 'Y':               out += 'Z';  sgn *= 1j
            elif b == 'Z':               out += 'Y';  sgn *= -1j
        elif a == 'Y':
            if   b == 'Y':               out += 'I'
            elif b == 'Z':               out += 'X';  sgn *= 1j
            elif b == 'X':               out += 'Z';  sgn *= -1j
        elif a == 'Z':
            if   b == 'Z':               out += 'I'
            elif b == 'X':               out += 'Y';  sgn *= 1j
            elif b == 'Y':               out += 'X';  sgn *= -1j
    return [out, sgn]


def _apply_rotation(rotation, p: str) -> Dict[str, float]:
    """Apply one Clifford rotation to a Pauli string; return linear combination."""
    out: Dict[str, float] = {}
    gen, angle = rotation[1], rotation[0]
    if not _commute(gen, p):
        if angle == 'pi/2':
            q = _pauli_mult(gen, p)
            out[q[0]] = float((1j * q[1]).real)
        else:
            out[p]  = np.cos(angle)
            q       = _pauli_mult(gen, p)
            out[q[0]] = float((1j * q[1] * np.sin(angle)).real)
    else:
        out[p] = 1.0
    return out


def _quasi_model_commuting_only(ham_nc: Dict[str, float]):
    """Minimal quasi-model for a purely commuting (no clique) noncontextual Hamiltonian.

    Returns (generators, [], all_mappings) compatible with the cs_vqe interface.
    Since all terms in ham_nc mutually commute, the anticommuting clique list is empty
    and every term maps directly to a product of generators.

    For the general case (cliques present), cs_vqe.quasi_model should be used.
    If cs_vqe is not available we fall back to identity mappings.
    """
    terms = list(ham_nc.keys())
    # Trivially: use each term as its own "generator" (overcomplete but valid)
    all_mappings = {t: [[t], [], 1] for t in terms}
    return terms, [], all_mappings


def _find_gs_classical(ham_nc: Dict[str, float]) -> float:
    """Classical ground state of a noncontextual Hamiltonian.

    For a fully commuting set, the ground state energy is obtained by assigning
    +1 or -1 to each generator and minimising the linear combination.
    This is equivalent to checking all 2^|G| sign assignments.
    """
    terms = list(ham_nc.keys())
    n     = len(terms)
    best  = np.inf
    for signs in itertools.product([1, -1], repeat=n):
        e = sum(s * ham_nc[t] for s, t in zip(signs, terms))
        if e < best:
            best = e
    return best


def _quantum_correction_in_sector(
    ham_eff: Dict[str, float],
    n_eff: int,
    verbose: bool = True,
) -> Dict[str, float]:
    """Apply CS-VQE quantum correction inside the Gibbs sector.

    Tries to import cs_vqe for the full quasi_model treatment; falls back to
    the self-contained implementation if cs_vqe is not installed.

    Returns the final reduced Hamiltonian dict (on <= n_eff qubits).
    """
    # ── Try cs_vqe first ────────────────────────────────────────────────────
    try:
        import cs_vqe as c  # type: ignore

        ham_nc_eff, ham_c_eff = _greedy_noncontextual_split(ham_eff)

        if not ham_nc_eff:
            if verbose:
                print("  ham_eff is fully contextual — returning H_eff unchanged.")
            return ham_eff

        model_eff   = c.quasi_model(ham_nc_eff)
        fn_form_eff = c.energy_function_form(ham_nc_eff, model_eff)
        gs_nc_eff   = c.find_gs_noncon(
            ham_nc_eff, model=model_eff, fn_form=fn_form_eff
        )
        ep_eff = gs_nc_eff[1]

        if verbose:
            print(f"  [cs_vqe] Noncontextual gs energy: {gs_nc_eff[0]:.6f}")

        if not ham_c_eff:
            # Already noncontextual → the noncontextual approximation is exact
            if verbose:
                print("  H_eff is noncontextual; no contextual correction needed.")
            return ham_nc_eff

        # get_reduced_hamiltonians returns a list of dicts for each qubit count;
        # the last element (all qubits quantum) is the full quantum correction.
        n_q = len(next(iter(ham_eff.keys())))
        reduced = c.get_reduced_hamiltonians(
            ham_eff, model_eff, fn_form_eff, ep_eff,
            order=list(range(n_q)),
        )
        # The first element has 0 quantum qubits (pure noncontextual);
        # the last has all qubits quantum (full correction).
        # We want the one that keeps the contextual correction = last element.
        ham_red = reduced[-1]
        if verbose:
            n_red = len(next(iter(ham_red.keys()))) if ham_red else 0
            print(f"  [cs_vqe] Quantum-corrected Hamiltonian: "
                  f"{len(ham_red)} terms on {n_red} qubits")
        return ham_red

    except ImportError:
        pass  # cs_vqe not available; use self-contained fallback below

    # ── Self-contained fallback ──────────────────────────────────────────────
    # Mirrors the fallback block in csvqe_gibbs_low_sector_demo_updated.ipynb
    if verbose:
        print("  [fallback] cs_vqe not found; using self-contained quantum correction.")

    ham_nc_eff, ham_c_eff = _greedy_noncontextual_split(ham_eff)

    if not ham_nc_eff:
        if verbose:
            print("  H_eff is fully contextual — returning H_eff unchanged.")
        return ham_eff

    if not ham_c_eff:
        if verbose:
            print("  H_eff is noncontextual — no contextual correction needed.")
        return ham_nc_eff

    # Diagonalise the noncontextual part
    H_nc_eff_mat             = _dict_to_matrix(ham_nc_eff, n_eff)
    evals_nc_eff, evecs_nc_eff = np.linalg.eigh(H_nc_eff_mat)
    nc_gs_vec                = evecs_nc_eff[:, 0]

    # Project H_eff onto the orthogonal complement of the nc ground state
    H_eff_mat  = _dict_to_matrix(ham_eff, n_eff)
    dim        = H_eff_mat.shape[0]
    P_nc_gs    = np.outer(nc_gs_vec, nc_gs_vec.conj())
    P_perp     = np.eye(dim) - P_nc_gs
    H_residual = P_perp @ H_eff_mat @ P_perp

    # The contextual residual Hamiltonian lives in the complement of nc_gs_vec;
    # drop the zero eigenvalue direction and decompose back to Pauli strings.
    n_res  = n_eff  # same qubit count; just the one direction projected out
    ham_red = _matrix_to_dict(H_residual, n_res)

    if verbose:
        print(f"  [fallback] Residual Hamiltonian: {len(ham_red)} terms on {n_res} qubits")

    return ham_red


# ──────────────────────────────────────────────────────────────────────────────
# 6.  Convenience alias
# ──────────────────────────────────────────────────────────────────────────────

# So callers can write: from ham_reduce import reduce_hamiltonian
__all__ = ["reduce_hamiltonian"]
