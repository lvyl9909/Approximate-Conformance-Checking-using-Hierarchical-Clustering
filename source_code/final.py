import xml.etree.ElementTree as ET
import pm4py
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
import numpy as np
import algorithm as clust_algorithm
from evaluation import eval_avg_leven
from pm4py.algo.conformance.alignments.petri_net import algorithm as alignments
from pm4py.objects.petri_net.importer import importer as pnml_importer
from pm4py.objects.log.exporter.xes import exporter as xes_exporter
from pm4py.objects.log.obj import EventLog, Trace, Event
from collections import defaultdict
import time

def parse_tracks_from_txt(file_path):
    """
    从指定的txt文件中解析中心迹所属的迹变体和出现次数最多的迹变体的轨迹代码集合。

    :param file_path: txt文件的路径
    :return: 包含中心迹所属的迹变体和出现次数最多的迹变体的轨迹代码列表
    """
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    center_tracks = set()
    most_frequent_tracks = set()

    is_center = False
    is_most_frequent = False

    for line in lines:
        line = line.strip()

        if line.startswith("簇") and "中心迹所属的迹变体索引" in line:
            is_center = True
            is_most_frequent = False
            continue

        if line.startswith("簇") and "出现次数最多的迹变体索引" in line:
            is_center = False
            is_most_frequent = True
            continue

        if line.startswith("迹变体") and "包含的轨迹" in line:
            continue

        if line:
            if is_center:
                center_tracks.add(line)
            elif is_most_frequent:
                most_frequent_tracks.add(line)

    return list(center_tracks), list(most_frequent_tracks)

def extract_traces_to_xml(xml_input_path, xml_output_center_path, xml_output_frequent_path, center_tracks_list,
                          most_frequent_tracks_list):
    """
    从指定的XML文件中提取指定的轨迹代码的trace节点，并生成两个新的XML文件。

    :param xml_input_path: 输入的XML文件路径
    :param xml_output_center_path: 输出的中心迹所属轨迹的XML文件路径
    :param xml_output_frequent_path: 输出的出现次数最多的轨迹的XML文件路径
    :param center_tracks_list: 中心迹所属的迹变体的轨迹代码列表
    :param most_frequent_tracks_list: 出现次数最多的迹变体的轨迹代码列表
    """
    tree = ET.parse(xml_input_path)
    root = tree.getroot()

    center_root = ET.Element(root.tag, root.attrib)
    most_frequent_root = ET.Element(root.tag, root.attrib)

    for child in root:
        if child.tag != 'trace':
            center_root.append(child)
            most_frequent_root.append(child)

    for trace in root.findall('trace'):
        trace_code = None
        for string in trace.findall('string'):
            if string.attrib.get('key') == 'concept:name':
                trace_code = string.attrib.get('value')
                break

        if trace_code in center_tracks_list:
            center_root.append(trace)
        if trace_code in most_frequent_tracks_list:
            most_frequent_root.append(trace)

    # 保存新的XML文件
    center_tree = ET.ElementTree(center_root)
    most_frequent_tree = ET.ElementTree(most_frequent_root)

    center_tree.write(xml_output_center_path, encoding='utf-8', xml_declaration=True)
    most_frequent_tree.write(xml_output_frequent_path, encoding='utf-8', xml_declaration=True)

    print(f"文件生成成功: {xml_output_center_path} 和 {xml_output_frequent_path}")

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

def execute_script():
    # Step 1: Generate cluster_out.txt
    with open("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2016/cluster_cp-10.txt", "w", encoding="utf-8") as f:
        log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2016/BPI2016_Questions.xes",
                             return_legacy_log_object=True)
        net, initial_marking, final_marking = pnml_importer.apply("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/petri_net-sepsis.pnml")

        trace_variants = generate_trace_variants(log)
        variant_id_map = {variant: idx + 1 for idx, variant in enumerate(trace_variants.keys())}

        y, variant_names = calculate_distance_matrix(trace_variants)
        dist_matrix = squareform(y)

        representative_traces_log = pm4py.objects.log.obj.EventLog(
            [trace_variants[variant][0] for variant in variant_names])

        start_time = time.time()
        cluster_tree, leafname, Z = clust_algorithm.apply(representative_traces_log,
                                                          "concept:name",
                                                          variant=clust_algorithm.Variants.VARIANT_AVG_LEVEN)

        num_clusters = 226
        clusters = fcluster(Z, num_clusters, criterion='maxclust')

        end_time = time.time()

        # Calculate the elapsed time in milliseconds
        clustering_time_ms = (end_time - start_time) * 1000
        print(f"层次聚类耗时: {clustering_time_ms:.3f} 毫秒")

        cluster_map = {i: [] for i in range(1, num_clusters + 1)}
        for variant_idx, cluster_id in enumerate(clusters):
            cluster_map[cluster_id].append(variant_names[variant_idx])

        for cluster_id, cluster_variants in cluster_map.items():
            medoid_variant, medoid_index = find_medoid_for_cluster(cluster_variants, dist_matrix, variant_names)
            f.write(f"簇 {cluster_id} 的中心迹所属的迹变体索引: {medoid_index + 1}\n")
            f.write(f"迹变体 {medoid_index + 1} 包含的轨迹:\n")
            for trace in trace_variants[medoid_variant]:
                f.write(f"  {trace.attributes['concept:name']}\n")
            f.write("\n")

            most_frequent_variant = find_most_frequent_variant(cluster_variants, trace_variants)
            most_frequent_variant_index = variant_names.index(most_frequent_variant)
            f.write(f"簇 {cluster_id} 中出现次数最多的迹变体索引: {most_frequent_variant_index + 1}\n")
            f.write(f"迹变体 {most_frequent_variant_index + 1} 包含的轨迹:\n")
            for trace in trace_variants[most_frequent_variant]:
                f.write(f"  {trace.attributes['concept:name']}\n")
            f.write("\n")

    # Step 2: Parse tracks from the generated txt file
    txt_file_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2016/cluster_cp-10.txt'

    center_tracks, most_frequent_tracks = parse_tracks_from_txt(txt_file_path)

    # Step 3: Extract traces from XML based on parsed track codes
    xml_input_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2016/BPI2016_Questions.xes'

    xml_output_center_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2016/incluster-medoid-10.xes'
    xml_output_frequent_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2016/incluster-frequency-10.xes'

    extract_traces_to_xml(
        xml_input_path,
        xml_output_center_path,
        xml_output_frequent_path,
        center_tracks,
        most_frequent_tracks
    )



if __name__ == "__main__":
    execute_script()
