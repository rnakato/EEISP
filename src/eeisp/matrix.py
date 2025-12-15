from __future__ import annotations

import multiprocessing as mp
from typing import Any, TypeAlias

import numpy as np
from numpy.typing import NDArray
from scipy.stats import binom

FloatArray: TypeAlias = NDArray[np.floating[Any]]
UInt16Array: TypeAlias = NDArray[np.uint16]
BoolArray: TypeAlias = NDArray[np.bool_]


def _argwrapper(args: tuple[Any, ...]) -> Any:
    func = args[0]
    fn = func if callable(func) else None
    if fn is None:
        raise TypeError("First arg must be callable")
    return fn(*args[1:])


def _calc_cdi_each_row(
    i: int,
    prob_array: FloatArray,
    count_array: NDArray[np.integer[Any]],
    ngene: int,
    ncell: int,
) -> FloatArray:
    _ = ngene
    x = count_array[i + 1 :]
    p = prob_array[i + 1 :]
    prob = binom.sf(x - 1, ncell, p)
    val = np.where(prob <= 0, -10_000_000.0, -(np.log10(prob)))
    return np.pad(val, (i + 1, 0), mode="constant", constant_values=0)


def _calc_eei_each_row(
    i: int,
    prob_array_i: FloatArray,
    count_array_i: NDArray[np.integer[Any]],
    prob_array_t_i: FloatArray,
    count_array_t_i: NDArray[np.integer[Any]],
    ngene: int,
    ncell: int,
) -> FloatArray:
    _ = ngene
    x1 = count_array_t_i[i + 1 :]
    p1 = prob_array_i[i + 1 :]
    prob1 = binom.sf(x1 - 1, ncell, p1)

    x2 = count_array_i[i + 1 :]
    p2 = prob_array_t_i[i + 1 :]
    prob2 = binom.sf(x2 - 1, ncell, p2)

    val = np.where(
        (prob1 <= 0) | (prob2 <= 0),
        -10_000_000.0,
        (-(np.log10(prob1)) + (-(np.log10(prob2)))) / 2,
    )
    return np.pad(val, (i + 1, 0), mode="constant", constant_values=0)


def _gen_matrix_multiprocess(
    prob_joint: FloatArray,
    count_joint: NDArray[np.integer[Any]],
    mat_type: str,
    ngene: int,
    ncell: int,
    *,
    ncore: int,
) -> FloatArray:
    ctx = mp.get_context("spawn")
    with ctx.Pool(ncore) as pool:
        func_args: list[tuple[Any, ...]] = []
        for i in range(0, ngene):
            if mat_type == "CDI":
                func_args.append((_calc_cdi_each_row, i, prob_joint[i], count_joint[i], ngene, ncell))
            elif mat_type == "EEI":
                func_args.append(
                    (
                        _calc_eei_each_row,
                        i,
                        prob_joint[i],
                        count_joint[i],
                        prob_joint.T[i],
                        count_joint.T[i],
                        ngene,
                        ncell,
                    )
                )
            else:
                raise ValueError(f"Illegal mat_type={mat_type!r}")

        results = pool.map(_argwrapper, func_args)

    matrix = np.array(results)
    matrix = matrix + matrix.T - np.diag(np.diag(matrix))
    return matrix


def _count_sum_nonzero_mat(i: int, is_nonzero_mat: BoolArray) -> NDArray[np.int64]:
    return np.sum(is_nonzero_mat[i, :] * is_nonzero_mat, axis=1)


def _count_sum_nonzero_mat_not_a(i: int, is_nonzero_mat: BoolArray, not_a: BoolArray) -> NDArray[np.int64]:
    return np.sum(is_nonzero_mat[i, :] * not_a, axis=1)


def generate_cdi_matrix(
    a: NDArray[np.number[Any]],
    *,
    threads: int,
    use_gpu: bool,
) -> FloatArray:
    ngene = int(a.shape[0])
    ncell = int(a.shape[1])

    is_nonzero_mat = a > 0
    p_nonzero = np.sum(is_nonzero_mat, axis=1) / ncell

    if use_gpu:
        import cupy as cp  # optional dependency; import only when requested

        p_nonzero_cp = cp.asarray(p_nonzero)
        is_nonzero_cp = cp.asarray(is_nonzero_mat)
        prob_joint_cp = cp.array(p_nonzero_cp * p_nonzero_cp[:, cp.newaxis], dtype="float32")
        count_joint_cp = cp.zeros((ngene, ngene), dtype="uint16")
        for i in range(ngene):
            count_joint_cp[i] = cp.sum(is_nonzero_cp[i] * is_nonzero_cp, axis=1)
        prob_joint = cp.asnumpy(prob_joint_cp)
        count_joint = cp.asnumpy(count_joint_cp)
    else:
        prob_joint = np.array(p_nonzero * p_nonzero[:, np.newaxis], dtype="float32")
        ctx = mp.get_context("spawn")
        with ctx.Pool(threads) as pool:
            func_args: list[tuple[Any, ...]] = []
            for i in range(ngene):
                func_args.append((_count_sum_nonzero_mat, i, is_nonzero_mat))
            count_joint_list = pool.map(_argwrapper, func_args)
        count_joint = np.array(count_joint_list, dtype=np.uint16)

    return _gen_matrix_multiprocess(prob_joint, count_joint, "CDI", ngene, ncell, ncore=threads)


def generate_eei_matrix(
    a: NDArray[np.number[Any]],
    *,
    threads: int,
    use_gpu: bool,
) -> FloatArray:
    ngene = int(a.shape[0])
    ncell = int(a.shape[1])

    is_nonzero_mat = a > 0
    p_nonzero = np.sum(is_nonzero_mat, axis=1) / ncell
    p_zero = np.sum(a == 0, axis=1) / ncell

    if use_gpu:
        import cupy as cp  # optional dependency; import only when requested

        p_nonzero_cp = cp.asarray(p_nonzero)
        p_zero_cp = cp.asarray(p_zero)
        is_nonzero_cp = cp.asarray(is_nonzero_mat)
        not_a_cp = cp.asarray(np.logical_not(a))
        prob_joint_cp = cp.array(p_nonzero_cp * p_zero_cp[:, cp.newaxis], dtype="float32")
        count_excl_cp = cp.zeros((ngene, ngene), dtype="uint16")
        for i in range(ngene):
            count_excl_cp[i] = cp.sum(is_nonzero_cp[i] * not_a_cp, axis=1)
        prob_joint = cp.asnumpy(prob_joint_cp)
        count_excl = cp.asnumpy(count_excl_cp)
    else:
        prob_joint = np.array(p_nonzero * p_zero[:, np.newaxis], dtype="float32")
        not_a = np.logical_not(a)
        ctx = mp.get_context("spawn")
        with ctx.Pool(threads) as pool:
            func_args: list[tuple[Any, ...]] = []
            for i in range(ngene):
                func_args.append((_count_sum_nonzero_mat_not_a, i, is_nonzero_mat, not_a))
            count_excl_list = pool.map(_argwrapper, func_args)
        count_excl = np.array(count_excl_list, dtype=np.uint16)

    return _gen_matrix_multiprocess(prob_joint, count_excl, "EEI", ngene, ncell, ncore=threads)


