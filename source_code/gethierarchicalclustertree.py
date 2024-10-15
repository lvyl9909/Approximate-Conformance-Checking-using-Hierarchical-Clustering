# Get Hierarchical Cluster Tree
import pm4py
from scipy.cluster.hierarchy import dendrogram, fcluster
import matplotlib.pyplot as plt
import algorithm as clust_algorithm
from merge_log import add_node, label_tree

def recursive_print_clusters(cluster_tree, rec_depth=0):
    if rec_depth > 0:
        print("\t" * rec_depth, cluster_tree["name"].split("-"))

    for child in cluster_tree["children"]:
        recursive_print_clusters(child, rec_depth + 1)

    if not cluster_tree["children"]:
        print("\t" * (rec_depth + 1), "END")

def recursive_print_clusters(cluster_tree, rec_depth=0):
    indent = "\t" * rec_depth
    node_name = cluster_tree["name"]
    distance = cluster_tree.get("distance", "N/A")

    print(f"{indent}{node_name} (distance: {distance})")

    for child in cluster_tree["children"]:
        recursive_print_clusters(child, rec_depth + 1)

    if not cluster_tree["children"]:
        print(f"{indent}  END")

def execute_script():
    log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis/sepsis.xes",
                         return_legacy_log_object=True)
    log = pm4py.sample_cases(log, num_cases=1050)

    # Perform hierarchical clustering on the 'concept:name' attribute of the log
    cluster_tree, Z, leafname = clust_algorithm.apply(log, "concept:name", variant=clust_algorithm.Variants.VARIANT_AVG_LEVEN)

    num_clusters = 10
    clusters = fcluster(Z, num_clusters, criterion='maxclust')

    # 将簇与对应的trace名称映射
    cluster_map = {}
    for trace, cluster_id in zip(log, clusters):
        trace_name = trace.attributes['concept:name']
        if cluster_id not in cluster_map:
            cluster_map[cluster_id] = []
        cluster_map[cluster_id].append(trace_name)

    # 打印簇及其对应的trace名称
    for cluster_id, trace_names in cluster_map.items():
        print(f"簇 {cluster_id}:")
        for trace_name in trace_names:
            print(f"  {trace_name}")

    # Print the cluster tree
    # recursive_print_clusters(cluster_tree)


if __name__ == "__main__":
    execute_script()