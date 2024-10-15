# import pm4py
# from scipy.cluster.hierarchy import fcluster
# from scipy.spatial.distance import squareform
# import numpy as np
# import algorithm as clust_algorithm
# from evaluation import eval_avg_leven
# import merge_log
# from collections import defaultdict
# from pm4py.algo.conformance.alignments.petri_net import algorithm as alignments
# from pm4py.objects.petri_net.importer import importer as pnml_importer
# from pm4py.objects.log.exporter.xes import exporter as xes_exporter
# from pm4py.objects.log.obj import EventLog, Trace, Event
# import os
#
#
# def generate_trace_variants(log):
#     """
#     Generate trace variants based on the activity sequence.
#     """
#     trace_variants = defaultdict(list)
#     for trace in log:
#         activity_sequence = "-".join(event['concept:name'] for event in trace)
#         trace_variants[activity_sequence].append(trace)
#     return trace_variants
#
#
# def calculate_distance_matrix(trace_variants):
#     """
#     Calculate the Levenshtein distance matrix for the given trace variants.
#     """
#     percent = 1
#     alpha = 0
#     variant_names = list(trace_variants.keys())
#     representative_traces = [trace_variants[variant][0] for variant in variant_names]
#     sublogs = [pm4py.objects.log.obj.EventLog([trace]) for trace in representative_traces]
#     y = eval_avg_leven(sublogs, percent, alpha)
#     return y, variant_names
#
#
# def find_medoid_for_cluster(cluster_variants, distance_matrix, variant_names):
#     """
#     Find the single medoid of a given cluster.
#     """
#     if len(cluster_variants) == 0:
#         return None
#     cluster_indices = [variant_names.index(variant) for variant in cluster_variants if variant in variant_names]
#     if not cluster_indices:
#         return None
#     sub_matrix = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
#     distance_sums = np.sum(sub_matrix, axis=1)
#     min_index = np.argmin(distance_sums)
#     medoid_index = cluster_indices[min_index]
#     medoid_variant = variant_names[medoid_index]
#     return medoid_variant, medoid_index
#
#
# def find_most_frequent_variant(cluster_variants, trace_variants):
#     """
#     Find the most frequent variant in a given cluster.
#     """
#     variant_frequency = {variant: len(trace_variants[variant]) for variant in cluster_variants}
#     most_frequent_variant = max(variant_frequency, key=variant_frequency.get)
#     return most_frequent_variant
#
#
# def align_trace_to_model(trace, net, initial_marking, final_marking):
#     """
#     Align a trace to the Petri net model to get the optimal alignment.
#     """
#     event_log = pm4py.objects.log.obj.EventLog([trace])
#     aligned_traces = alignments.apply(event_log, net, initial_marking, final_marking)
#     return aligned_traces
#
#
# def execute_script():
#     # 检查并创建输出目录
#     output_dir = "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis"
#     if not os.path.exists(output_dir):
#         os.makedirs(output_dir)
#
#     # 打开文件以写入输出
#     with open(os.path.join(output_dir, "test.txt"), "w", encoding="utf-8") as f:
#         # 读取日志和PNML文件
#         log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis/sepsis.xes",
#                              return_legacy_log_object=True)
#         net, initial_marking, final_marking = pnml_importer.apply("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis/spesis_0.4.pnml")
#
#         # 生成迹变体
#         trace_variants = generate_trace_variants(log)
#         variant_id_map = {variant: idx + 1 for idx, variant in enumerate(trace_variants.keys())}
#
#         # 计算迹变体之间的Levenshtein距离矩阵，只使用代表性的轨迹
#         y, variant_names = calculate_distance_matrix(trace_variants)
#         dist_matrix = squareform(y)
#
#         # 对迹变体进行层次聚类，只使用代表性的轨迹
#         representative_traces_log = pm4py.objects.log.obj.EventLog(
#             [trace_variants[variant][0] for variant in variant_names])
#         cluster_tree, leafname, Z = clust_algorithm.apply(representative_traces_log,
#                                                           "concept:name",
#                                                           variant=clust_algorithm.Variants.VARIANT_AVG_LEVEN)
#
#         num_clusters = 456
#         clusters = fcluster(Z, num_clusters, criterion='maxclust')
#
#         # 将簇与对应的迹变体名称映射
#         cluster_map = {i: [] for i in range(1, num_clusters + 1)}
#         for variant_idx, cluster_id in enumerate(clusters):
#             cluster_map[cluster_id].append(variant_names[variant_idx])
#
#         # 遍历每个簇并找到medoid，并生成最优对齐
#         for cluster_id, cluster_variants in cluster_map.items():
#             # 添加检查，避免 None 解包错误
#             medoid_info = find_medoid_for_cluster(cluster_variants, dist_matrix, variant_names)
#             if medoid_info is None:
#                 f.write(f"簇 {cluster_id} 没有可用的中心迹变体。\n")
#                 continue
#
#             medoid_variant, medoid_index = medoid_info
#             f.write(f"簇 {cluster_id} 的中心迹所属的迹变体索引: {medoid_index + 1}\n")
#             f.write(f"原迹变体 {medoid_index + 1} 包含的轨迹:\n")
#             for trace in trace_variants[medoid_variant]:
#                 f.write(f"  {trace.attributes['concept:name']}\n")
#
#             # Align the medoid trace to the model
#             medoid_trace = trace_variants[medoid_variant][0]
#             aligned_traces = align_trace_to_model(medoid_trace, net, initial_marking, final_marking)
#             f.write(f"对齐后的中心迹:\n")
#             for aligned_trace in aligned_traces:
#                 for move in aligned_trace['alignment']:
#                     if move[1] and move[1] != "None":  # Only print the event name and skip 'None'
#                         f.write(f"  {move[1]}\n")
#             f.write("\n")
#
#             # 查找并打印簇中出现次数最多的迹变体
#             most_frequent_variant = find_most_frequent_variant(cluster_variants, trace_variants)
#             most_frequent_variant_index = variant_names.index(most_frequent_variant)
#             f.write(f"簇 {cluster_id} 中出现次数最多的迹变体索引: {most_frequent_variant_index + 1}\n")
#             f.write(f"原迹变体 {most_frequent_variant_index + 1} 包含的轨迹:\n")
#             for trace in trace_variants[most_frequent_variant]:
#                 f.write(f"  {trace.attributes['concept:name']}\n")
#
#             # Align the most frequent trace to the model
#             most_frequent_trace = trace_variants[most_frequent_variant][0]
#             aligned_traces_frequent = align_trace_to_model(most_frequent_trace, net, initial_marking, final_marking)
#             f.write(f"对齐后的最频繁迹:\n")
#             for aligned_trace in aligned_traces_frequent:
#                 for move in aligned_trace['alignment']:
#                     if move[1] and move[1] != "None":  # Only print the event name and skip 'None'
#                         f.write(f"  {move[1]}\n")
#             f.write("\n")
#
#
# if __name__ == "__main__":
#     execute_script()

import pm4py
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
import numpy as np
import algorithm as clust_algorithm
from evaluation import eval_avg_leven
import merge_log
from collections import defaultdict
from pm4py.algo.conformance.alignments.petri_net import algorithm as alignments
from pm4py.objects.petri_net.importer import importer as pnml_importer
from pm4py.objects.log.exporter.xes import exporter as xes_exporter
from pm4py.objects.log.obj import EventLog, Trace, Event
import os


def generate_trace_variants(log):
    """
    Generate trace variants based on the activity sequence.
    """
    trace_variants = defaultdict(list)
    for trace in log:
        activity_sequence = "-".join(event['concept:name'] for event in trace)
        trace_variants[activity_sequence].append(trace)
    return trace_variants


def calculate_distance_matrix(trace_variants):
    """
    Calculate the Levenshtein distance matrix for the given trace variants.
    """
    percent = 1
    alpha = 0
    variant_names = list(trace_variants.keys())
    representative_traces = [trace_variants[variant][0] for variant in variant_names]
    sublogs = [pm4py.objects.log.obj.EventLog([trace]) for trace in representative_traces]
    y = eval_avg_leven(sublogs, percent, alpha)
    return y, variant_names


def find_medoid_for_cluster(cluster_variants, distance_matrix, variant_names):
    """
    Find the single medoid of a given cluster.
    """
    if len(cluster_variants) == 0:
        return None
    cluster_indices = [variant_names.index(variant) for variant in cluster_variants if variant in variant_names]
    if not cluster_indices:
        return None
    sub_matrix = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
    distance_sums = np.sum(sub_matrix, axis=1)
    min_index = np.argmin(distance_sums)
    medoid_index = cluster_indices[min_index]
    medoid_variant = variant_names[medoid_index]
    return medoid_variant, medoid_index

def find_most_frequent_variant(cluster_variants, trace_variants):
    """
    Find the most frequent variant in a given cluster.
    """
    variant_frequency = {variant: len(trace_variants[variant]) for variant in cluster_variants}
    most_frequent_variant = max(variant_frequency, key=variant_frequency.get)
    return most_frequent_variant

def align_trace_to_model(trace, net, initial_marking, final_marking):
    """
    Align a trace to the Petri net model to get the optimal alignment.
    """
    event_log = pm4py.objects.log.obj.EventLog([trace])
    aligned_traces = alignments.apply(event_log, net, initial_marking, final_marking)
    return aligned_traces

def execute_script():
    # 检查并创建输出目录
    output_dir = "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 打开文件以写入输出
    with open(os.path.join(output_dir, "test.txt"), "w", encoding="utf-8") as f:
        # 读取日志和PNML文件
        log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis/sepsis.xes",
                             return_legacy_log_object=True)
        net, initial_marking, final_marking = pnml_importer.apply("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis/spesis_0.4.pnml")

        # 生成迹变体
        trace_variants = generate_trace_variants(log)
        variant_id_map = {variant: idx + 1 for idx, variant in enumerate(trace_variants.keys())}

        # 计算迹变体之间的Levenshtein距离矩阵，只使用代表性的轨迹
        y, variant_names = calculate_distance_matrix(trace_variants)
        dist_matrix = squareform(y)

        # 对迹变体进行层次聚类，只使用代表性的轨迹
        representative_traces_log = pm4py.objects.log.obj.EventLog(
            [trace_variants[variant][0] for variant in variant_names])
        cluster_tree, leafname, Z = clust_algorithm.apply(representative_traces_log,
                                                          "concept:name",
                                                          variant=clust_algorithm.Variants.VARIANT_AVG_LEVEN)

        num_clusters = 456
        clusters = fcluster(Z, num_clusters, criterion='maxclust')

        # 将簇与对应的迹变体名称映射
        cluster_map = {i: [] for i in range(1, num_clusters + 1)}
        for variant_idx, cluster_id in enumerate(clusters):
            cluster_map[cluster_id].append(variant_names[variant_idx])

        # 遍历每个簇并找到medoid，并生成最优对齐
        for cluster_id, cluster_variants in cluster_map.items():
            # 添加检查，避免 None 解包错误
            medoid_info = find_medoid_for_cluster(cluster_variants, dist_matrix, variant_names)
            if medoid_info is None:
                f.write(f"簇 {cluster_id} 没有可用的中心迹变体。\n")
                continue

            medoid_variant, medoid_index = medoid_info
            f.write(f"簇 {cluster_id} 的中心迹所属的迹变体索引: {medoid_index + 1}\n")
            f.write(f"原迹变体 {medoid_index + 1} 包含的轨迹:\n")
            for trace in trace_variants[medoid_variant]:
                f.write(f"  {trace.attributes['concept:name']}\n")

            # Align the medoid trace to the model
            medoid_trace = trace_variants[medoid_variant][0]
            aligned_traces = align_trace_to_model(medoid_trace, net, initial_marking, final_marking)
            f.write(f"对齐后的中心迹:\n")
            for aligned_trace in aligned_traces:
                for move in aligned_trace['alignment']:
                    # Print both the symbol (">>") and the event name
                    f.write(f"  {move[0]} {move[1]}\n")
            f.write("\n")

            # 查找并打印簇中出现次数最多的迹变体
            most_frequent_variant = find_most_frequent_variant(cluster_variants, trace_variants)
            most_frequent_variant_index = variant_names.index(most_frequent_variant)
            f.write(f"簇 {cluster_id} 中出现次数最多的迹变体索引: {most_frequent_variant_index + 1}\n")
            f.write(f"原迹变体 {most_frequent_variant_index + 1} 包含的轨迹:\n")
            for trace in trace_variants[most_frequent_variant]:
                f.write(f"  {trace.attributes['concept:name']}\n")

            # Align the most frequent trace to the model
            most_frequent_trace = trace_variants[most_frequent_variant][0]
            aligned_traces_frequent = align_trace_to_model(most_frequent_trace, net, initial_marking, final_marking)
            f.write(f"对齐后的最频繁迹:\n")
            for aligned_trace in aligned_traces_frequent:
                for move in aligned_trace['alignment']:
                    # Print both the symbol (">>") and the event name
                    f.write(f"  {move[0]} {move[1]}\n")
            f.write("\n")

if __name__ == "__main__":
    execute_script()