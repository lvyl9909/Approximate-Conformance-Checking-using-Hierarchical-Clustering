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
    # 打开文件以写入输出
    with open("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2012/cluster_out.txt", "w", encoding="utf-8") as f:
        # 读取日志和PNML文件
        log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2012/BPIC2012.xes",
                             return_legacy_log_object=True)
        net, initial_marking, final_marking = pnml_importer.apply("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/petri_net-sepsis.pnml")

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
            medoid_variant, medoid_index = find_medoid_for_cluster(cluster_variants, dist_matrix, variant_names)
            f.write(f"簇 {cluster_id} 的中心迹所属的迹变体索引: {medoid_index + 1}\n")
            f.write(f"迹变体 {medoid_index + 1} 包含的轨迹:\n")
            for trace in trace_variants[medoid_variant]:
                f.write(f"  {trace.attributes['concept:name']}\n")
            f.write("\n")

            # 查找并打印簇中出现次数最多的迹变体
            most_frequent_variant = find_most_frequent_variant(cluster_variants, trace_variants)
            most_frequent_variant_index = variant_names.index(most_frequent_variant)
            f.write(f"簇 {cluster_id} 中出现次数最多的迹变体索引: {most_frequent_variant_index + 1}\n")
            f.write(f"迹变体 {most_frequent_variant_index + 1} 包含的轨迹:\n")
            for trace in trace_variants[most_frequent_variant]:
                f.write(f"  {trace.attributes['concept:name']}\n")
            f.write("\n")

if __name__ == "__main__":
    execute_script()
#
# def execute_script():
#     # log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/Lfull_final.xes",
#     #                      return_legacy_log_object=True)
#     log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis/sepsis.xes",
#                          return_legacy_log_object=True)
#     # 读取PNML文件，加载Petri网模型
#     # net, initial_marking, final_marking = pnml_importer.apply(
#     #     "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/ourpetrinet.pnml")
#     net, initial_marking, final_marking = pnml_importer.apply("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/petri_net-sepsis.pnml")
#     # 生成迹变体
#     trace_variants = generate_trace_variants(log)
#     variant_id_map = {variant: idx + 1 for idx, variant in enumerate(trace_variants.keys())}
#
#     # 打印每个迹变体编号和内部的trace concept:name
#     for variant, traces in trace_variants.items():
#         print(f"Trace Variant {variant_id_map[variant]} (Activity Sequence: {variant}):")
#         for trace in traces:
#             print(f"  Trace concept:name: {trace.attributes['concept:name']}")
#         print()
#     # 计算迹变体之间的Levenshte因距离矩阵，只使用代表性的轨迹
#     y, variant_names = calculate_distance_matrix(trace_variants)
#     dist_matrix = squareform(y)
#
#     # 对迹变体进行层次聚类，只使用代表性的轨迹
#     representative_traces_log = pm4py.objects.log.obj.EventLog(
#         [trace_variants[variant][0] for variant in variant_names])
#     cluster_tree, leafname, Z = clust_algorithm.apply(representative_traces_log,
#                                                       "concept:name",
#                                                       variant=clust_algorithm.Variants.VARIANT_AVG_LEVEN)
#
#     num_clusters = 105
#     clusters = fcluster(Z, num_clusters, criterion='maxclust')
#
#     # 将簇与对应的迹变体名称映射
#     cluster_map = {i: [] for i in range(1, num_clusters + 1)}
#     for variant_idx, cluster_id in enumerate(clusters):
#         cluster_map[cluster_id].append(variant_names[variant_idx])
#     for cluster_id, variant_names_in_cluster in cluster_map.items():
#         # print(f"簇 {cluster_id}:")
#         for variant_name in variant_names_in_cluster:
#             trace_count = len(trace_variants[variant_name])
#             # print(f"  Trace Variant: {variant_name} (包含 {trace_count} 条轨迹)")
#     # 创建一个新的事件日志用于保存对齐后的轨迹
#     aligned_log = EventLog()
#     most_frequent_aligned_log = EventLog()
#
#     # 遍历每个簇并找到medoid，并生成最优对齐
#     for cluster_id, cluster_variants in cluster_map.items():
#         medoid_variant, medoid_index = find_medoid_for_cluster(cluster_variants, dist_matrix, variant_names)
#         print(f"簇 {cluster_id} 的中心迹所属的迹变体索引: {medoid_index + 1}")
#         print(f"迹变体 {medoid_index + 1} 包含的轨迹:")
#         for trace in trace_variants[medoid_variant]:
#             print(f"  {trace.attributes['concept:name']}")
#         print()
#
#         # # 使用对齐算法将中心迹对齐到模型
#         # trace_to_align = trace_variants[medoid_variant][0]
#         # aligned_traces = align_trace_to_model(trace_to_align, net, initial_marking, final_marking)
#         #
#         # # 打印对齐结果并将其添加到对齐日志中（保持原样）
#         # for aligned_trace in aligned_traces:
#         #     aligned_trace_log = Trace()
#         #     print(f"对齐结果 (Trace):")
#         #     for move in aligned_trace['alignment']:
#         #         if move[0] == '>>':
#         #             # 模型移动，记录并打印
#         #             event = Event({"concept:name": move[1]})
#         #             print(f"  模型移动: {move[1]}")
#         #         else:
#         #             # 日志同步或日志移动，记录并打印
#         #             event = Event({"concept:name": move[1]})
#         #             print(f"  日志事件: {move[1]}")
#         #         aligned_trace_log.append(event)
#         #     aligned_log.append(aligned_trace_log)
#
#         # 查找并打印簇中出现次数最多的迹变体
#         most_frequent_variant = find_most_frequent_variant(cluster_variants, trace_variants)
#         most_frequent_variant_index = variant_names.index(most_frequent_variant)
#         print(f"簇 {cluster_id} 中出现次数最多的迹变体索引: {most_frequent_variant_index + 1}")
#         print(f"迹变体 {most_frequent_variant_index + 1} 包含的轨迹:")
#         for trace in trace_variants[most_frequent_variant]:
#             print(f"  {trace.attributes['concept:name']}")
#         print()
#     #
#     #     # 使用对齐算法将最频繁迹变体对齐到模型
#     #     trace_to_align = trace_variants[most_frequent_variant][0]
#     #     aligned_traces_frequent = align_trace_to_model(trace_to_align, net, initial_marking, final_marking)
#     #
#     #     for aligned_trace in aligned_traces_frequent:
#     #         aligned_trace_log = Trace()
#     #         print(f"对齐结果 (Trace):")
#     #         for move in aligned_trace['alignment']:
#     #             if move[0] == '>>':
#     #                 # 模型移动，记录并打印
#     #                 event = Event({"concept:name": move[1]})
#     #                 print(f"  模型移动: {move[1]}")
#     #             else:
#     #                 # 日志同步或日志移动，记录并打印
#     #                 event = Event({"concept:name": move[1]})
#     #                 print(f"  日志事件: {move[1]}")
#     #             aligned_trace_log.append(event)
#     #         most_frequent_aligned_log.append(aligned_trace_log)
#     # # 创建一个新的事件日志用于导出，仅保留实际活动
#     # filtered_log = EventLog()
#     # for trace in aligned_log:
#     #     filtered_trace = Trace()
#     #     for event in trace:
#     #         # 确保事件对象具有 "concept:name" 属性且其值不是 None
#     #         if event["concept:name"] is not None and event["concept:name"].isalpha():
#     #             filtered_trace.append(event)
#     #     filtered_log.append(filtered_trace)
#     #
#     # filtered_log_frequent = EventLog()
#     # for trace in most_frequent_aligned_log:
#     #     filtered_trace = Trace()
#     #     for event in trace:
#     #         # 确保事件对象具有 "concept:name" 属性且其值不是 None
#     #         if event["concept:name"] is not None and event["concept:name"].isalpha():
#     #             filtered_trace.append(event)
#     #     filtered_log_frequent.append(filtered_trace)
#     #
#     # # 将对齐后的实际活动日志导出为XES文件
#     # xes_exporter.apply(filtered_log, "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/aligned_traces_final.xes")
#     # xes_exporter.apply(filtered_log_frequent, "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/aligned_traces_frequent_final.xes")
#
#
# if __name__ == "__main__":
#     execute_script()

