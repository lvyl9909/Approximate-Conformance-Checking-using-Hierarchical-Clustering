'''
    This file is part of PM4Py (More Info: https://pm4py.fit.fraunhofer.de).

    PM4Py is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PM4Py is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PM4Py.  If not, see <https://www.gnu.org/licenses/>.
'''
from scipy.spatial.distance import squareform
import numpy as np
from pm4py.algo.clustering.trace_attribute_driven.variants import act_dist_calc
from pm4py.algo.clustering.trace_attribute_driven.variants import suc_dist_calc
import leven_dist_calc
from pm4py.algo.clustering.trace_attribute_driven.dfg import dfg_dist


def dfg_dis(loglist, percent, alpha):
    size = len(loglist)
    dist_mat = np.zeros((size, size))

    for i in range(0, size - 1):
        for j in range(i + 1, size):
            (dist_act, dist_dfg) = dfg_dist.dfg_dist_calc(loglist[i], loglist[j])
            dist_mat[i][j] = dist_act * alpha + dist_dfg * (1 - alpha)
            dist_mat[j][i] = dist_mat[i][j]
    y = squareform(dist_mat)
    return y


def eval_avg_variant(loglist, percent, alpha):
    size = len(loglist)
    dist_mat = np.zeros((size, size))

    for i in range(0, size - 1):
        for j in range(i + 1, size):
            dist_act = act_dist_calc.act_sim_percent_avg(loglist[i], loglist[j], percent, percent)
            dist_suc = suc_dist_calc.suc_sim_percent_avg(loglist[i], loglist[j], percent, percent)
            dist_mat[i][j] = dist_act * alpha + dist_suc * (1 - alpha)
            dist_mat[j][i] = dist_mat[i][j]
    y = squareform(dist_mat)

    return y


def eval_DMM_variant(loglist, percent, alpha):
    size = len(loglist)
    dist_mat = np.zeros((size, size))

    for i in range(0, size - 1):
        for j in range(i + 1, size):
            dist_act = act_dist_calc.act_sim_percent(loglist[i], loglist[j], percent, percent)
            dist_suc = suc_dist_calc.suc_sim_percent(loglist[i], loglist[j], percent, percent)
            dist_mat[i][j] = dist_act * alpha + dist_suc * (1 - alpha)
            dist_mat[j][i] = dist_mat[i][j]
    y = squareform(dist_mat)
    return y
    # for i in range(0, size - 1):
    #     for j in range(i + 1, size):
    #         try:
    #             dist_mat[i][j] = leven_dist_calc.leven_dist_avg(loglist[i], loglist[j], percent, percent)
    #             dist_mat[j][i] = dist_mat[i][j]
    #         except IndexError as e:
    #             print(f"IndexError: {e} at i={i}, j={j}")
    #             print(f"loglist[i]: {loglist[i]}")
    #             print(f"loglist[j]: {loglist[j]}")
    #             dist_mat[i][j] = np.inf  # Assign a large distance value to handle the error
    #             dist_mat[j][i] = dist_mat[i][j]
    # y = squareform(dist_mat)
    # return y

def eval_avg_leven(loglist, percent, alpha):
    size = len(loglist)
    dist_mat = np.zeros((size, size))

    for i in range(0, size - 1):
        for j in range(i + 1, size):
            dist_mat[i][j] = leven_dist_calc.leven_dist_avg(loglist[i], loglist[j], percent, percent)
            dist_mat[j][i] = dist_mat[i][j]
    y = squareform(dist_mat)
    return y


def eval_DMM_leven(loglist, percent, alpha):
    size = len(loglist)
    dist_mat = np.zeros((size, size))

    for i in range(0, size - 1):
        for j in range(i + 1, size):
            dist_mat[i][j] = leven_dist_calc.leven_dist(loglist[i], loglist[j], percent, percent)
            dist_mat[j][i] = dist_mat[i][j]
    y = squareform(dist_mat)
    return y

def eval_upgma_variant(loglist, percent, alpha):
    size = len(loglist)
    dist_mat = np.zeros((size, size))

    for i in range(0, size - 1):
        for j in range(i + 1, size):
            dist_act = act_dist_calc.act_sim_percent_avg(loglist[i], loglist[j], percent, percent)
            dist_suc = suc_dist_calc.suc_sim_percent_avg(loglist[i], loglist[j], percent, percent)
            dist_mat[i][j] = dist_act * alpha + dist_suc * (1 - alpha)
            dist_mat[j][i] = dist_mat[i][j]

    y = squareform(dist_mat)
    return y

def eval_upgma_leven(loglist, percent, alpha):
    size = len(loglist)
    dist_mat = np.zeros((size, size))

    for i in range(0, size - 1):
        for j in range(i + 1, size):
            dist_mat[i][j] = leven_dist_calc.leven_dist_avg(loglist[i], loglist[j], percent, percent)
            dist_mat[j][i] = dist_mat[i][j]

    y = squareform(dist_mat)
    return y

def eval_raw_leven(loglist, percent, alpha):
    size = len(loglist)
    dist_mat = np.zeros((size, size))

    for i in range(0, size - 1):
        for j in range(i + 1, size):
            dist_mat[i][j] = leven_dist_calc.leven_dist_raw(loglist[i], loglist[j], percent, percent)
            dist_mat[j][i] = dist_mat[i][j]
    y = squareform(dist_mat)
    return y

def calculate_frequency(trace_variants):
    """
    计算每个轨迹变体的出现频率。

    :param trace_variants: 包含每个轨迹变体及其对应轨迹的字典，键是轨迹变体的名称，值是轨迹的列表。
    :return: 包含每个轨迹变体出现次数的字典。
    """
    frequency_dict = {}

    # 遍历所有轨迹变体，计算每个轨迹变体的频率
    for variant, traces in trace_variants.items():
        frequency_dict[variant] = len(traces)  # 轨迹的数量即为该变体的频率

    return frequency_dict

def weighted_distance(loglist, percent):
    size = len(loglist)
    dist_mat = np.zeros((size, size))

    for i in range(0, size - 1):
        for j in range(i + 1, size):
            f_v1 = calculate_frequency(loglist[i])  # 假设有一个函数计算轨迹的频率
            f_v2 = calculate_frequency(loglist[j])  # 同上
            d_N = leven_dist_calc.leven_dist(loglist[i], loglist[j], percent, percent)  # 利用已有的Levenshtein距离计算

            # 根据给定的公式计算加权距离
            dist_weighted = (f_v1 * f_v2 * d_N) / max(f_v1 ** 2, f_v2 ** 2)
            dist_mat[i][j] = dist_weighted
            dist_mat[j][i] = dist_weighted

    y = squareform(dist_mat)
    return y
