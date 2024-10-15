
# import pm4py
# from scipy.cluster.hierarchy import fcluster
# from scipy.spatial.distance import squareform
# import numpy as np
# import algorithm as clust_algorithm
# from evaluation import eval_avg_leven
# import merge_log
#
# def calculate_distance_matrix(log):
#     """
#     Calculate the Levenshtein distance matrix for the given traces.
#     """
#     percent = 1
#     alpha = 0
#     list_of_vals = [trace.attributes["concept:name"] for trace in log]
#
#     list_log = [merge_log.log2sublog(log, val, "concept:name") for val in list_of_vals]
#
#     y = eval_avg_leven(list_log, percent, alpha)
#
#     return y, list_of_vals
#
# def find_medoid_for_cluster(cluster_traces, distance_matrix, list_of_vals):
#     """
#     Find the medoid(s) of a given cluster.
#     """
#     if len(cluster_traces) == 0:
#         return None
#
#     # Get indices of the traces in the cluster
#     cluster_indices = [list_of_vals.index(trace) for trace in cluster_traces]
#
#     # Extract sub-matrix for the current cluster
#     sub_matrix = distance_matrix[np.ix_(cluster_indices, cluster_indices)]
#
#     # Calculate sum of distances for each trace in the cluster
#     distance_sums = np.sum(sub_matrix, axis=1)
#
#     # Find the index with the minimum distance sum
#     min_index = np.argmin(distance_sums)
#
#     # Optionally, find the second closest trace if needed
#     medoid_indices = [cluster_indices[min_index]]
#
#     if len(cluster_indices) > 1:
#         second_min_index = np.argsort(distance_sums)[1]
#         medoid_indices.append(cluster_indices[second_min_index])
#
#     # Get the trace names for the medoid(s)
#     medoids = [list_of_vals[idx] for idx in medoid_indices]
#     return medoids
#
# def execute_script():
#     log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/Lfull_final.xes",
#                          return_legacy_log_object=True)
#     log = pm4py.sample_cases(log, num_cases=1050)
#
#     # 计算Levenshtein距离矩阵
#     y, list_of_vals = calculate_distance_matrix(log)
#
#     # 将压缩的距离矩阵转换为方阵
#     dist_matrix = squareform(y)
#
#     # 对日志进行层次聚类
#     cluster_tree, leafname, Z = clust_algorithm.apply(log, "concept:name", variant=clust_algorithm.Variants.VARIANT_AVG_LEVEN)
#
#     num_clusters = 3
#     clusters = fcluster(Z, num_clusters, criterion='maxclust')
#
#     # 将簇与对应的trace名称映射
#     cluster_map = {i: [] for i in range(1, num_clusters + 1)}
#     for trace, cluster_id in zip(log, clusters):
#         trace_name = trace.attributes['concept:name']
#         cluster_map[cluster_id].append(trace_name)
#
#     # 打印簇及其对应的trace名称
#     for cluster_id, trace_names in cluster_map.items():
#         print(f"簇 {cluster_id}:")
#         for trace_name in trace_names:
#             print(f"  {trace_name}")
#
#     # 遍历每个簇并找到medoid(s)
#     for cluster_id, cluster_traces in cluster_map.items():
#         medoids = find_medoid_for_cluster(cluster_traces, dist_matrix, list_of_vals)
#         print(f"簇 {cluster_id} 的 Medoid(s): {medoids}")
#
# if __name__ == "__main__":
#     execute_script()

# def recursive_print_clusters(cluster_tree, rec_depth=0):
#     if rec_depth > 0:
#         print("\t" * rec_depth, cluster_tree["name"].split("-"))
#
#     for child in cluster_tree["children"]:
#         recursive_print_clusters(child, rec_depth + 1)
#
#     if not cluster_tree["children"]:
#         print("\t" * (rec_depth + 1), "END")
#
# def recursive_print_clusters(cluster_tree, rec_depth=0):
#     indent = "\t" * rec_depth
#     node_name = cluster_tree["name"]
#     distance = cluster_tree.get("distance", "N/A")
#
#     print(f"{indent}{node_name} (distance: {distance})")
#
#     for child in cluster_tree["children"]:
#         recursive_print_clusters(child, rec_depth + 1)
#
#     if not cluster_tree["children"]:
#         print(f"{indent}  END")
#
# def execute_script():
#     log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis/sepsis.xes",
#                          return_legacy_log_object=True)
#     log = pm4py.sample_cases(log, num_cases=1050)
#
#     # Perform hierarchical clustering on the 'concept:name' attribute of the log
#     cluster_tree, leafname, Z = clust_algorithm.apply(log, "concept:name", variant=clust_algorithm.Variants.VARIANT_AVG_LEVEN)
#
#     num_clusters = 10
#     clusters = fcluster(Z, num_clusters, criterion='maxclust')
#
#     # 将簇与对应的trace名称映射
#     cluster_map = {}
#     for trace, cluster_id in zip(log, clusters):
#         trace_name = trace.attributes['concept:name']
#         if cluster_id not in cluster_map:
#             cluster_map[cluster_id] = []
#         cluster_map[cluster_id].append(trace_name)
#
#     # 打印簇及其对应的trace名称
#     for cluster_id, trace_names in cluster_map.items():
#         print(f"簇 {cluster_id}:")
#         for trace_name in trace_names:
#             print(f"  {trace_name}")
#
#     # Iterate through each cluster and find medoid(s)
#     for cluster_id, cluster_traces in cluster_map.items():
#         medoids = find_medoid_for_cluster(cluster_traces, dist_matrix, list_of_vals)
#         print(f"Cluster {cluster_id} Medoid(s): {medoids}")
#     # Print the cluster tree
#     # recursive_print_clusters(cluster_tree)
#
#
# if __name__ == "__main__":
#     execute_script()










