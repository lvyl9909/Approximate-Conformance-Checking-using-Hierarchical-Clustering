
# K-means method
# K-center method
# Inter-cluster medoid method
# Inter-cluster high-frequency method
# Inter-cluster multialignment method
# Inter-cluster consensus path method

import pm4py
import merge_log
from evaluation import eval_avg_leven
import numpy as np
from pm4py.objects.log.obj import EventLog
from scipy.spatial.distance import squareform
from pm4py.objects.conversion.log import converter as log_converter
from collections import Counter
import random
from pm4py.objects.log.importer.xes import importer as xes_importer
from pm4py.objects.petri_net.importer import importer as petri_importer
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
from discounted_a_star import apply as multii, Parameters
from difflib import SequenceMatcher
from scipy.spatial.distance import cdist

def filter_logs(log, keys):
    """
    Filter out cases where concept:name is in keys, retaining only one trace per case.
    """
    filtered_logs = []
    case_ids = set()
    for trace in log:
        if trace.attributes["concept:name"] in keys and trace.attributes["concept:name"] not in case_ids:
            filtered_logs.append(trace)
            case_ids.add(trace.attributes["concept:name"])
    return EventLog(filtered_logs)


def calculate_distance_matrix(log):
    """
    Calculate the Levenshtein distance matrix for the given traces.
    """
    percent = 1
    alpha = 0
    list_of_vals = [trace.attributes["concept:name"] for trace in log]

    list_log = [merge_log.log2sublog(log, val, "concept:name") for val in list_of_vals]

    y = eval_avg_leven(list_log, percent, alpha)

    print(f"Generated Distance Matrix y: {y}")  # Print y to see its actual content

    return y, list_of_vals


def find_case_with_highest_frequency(log, keys):
    """
    Find the case with the highest number of traces in the given log.
    log: The event log.
    keys: The set of case IDs to consider.
    Returns the case ID with the highest frequency of traces.
    """
    max_traces = 0
    case_with_max_traces = None

    for trace in log:
        case_id = trace.attributes["concept:name"]
        if case_id in keys:
            trace_count = len(log_converter.apply(trace))
            if trace_count > max_traces:
                max_traces = trace_count
                case_with_max_traces = case_id

    return case_with_max_traces


def find_new_medoid(y, list_of_vals):
    """
    Select the new medoid based on the Levenshtein distance matrix.
    y: The condensed one-dimensional distance matrix.
    list_of_vals: The list of trace names.
    """
    # Convert the condensed one-dimensional distance matrix back to a square form
    dist_matrix = squareform(y)

    # Calculate the sum of distances from each trace to all other traces
    dist_sums = np.sum(dist_matrix, axis=1)

    # Find the index of the trace with the smallest sum of distances
    min_index = np.argmin(dist_sums)

    # Return the corresponding trace name
    return list_of_vals[min_index]

##K-means
def calculate_edit_distance(trace1, trace2):
    # 使用difflib中的SequenceMatcher计算编辑距离
    return len(trace1) + len(trace2) - 2 * sum(triple.size for triple in SequenceMatcher(None, trace1, trace2).get_matching_blocks())

def get_trace_representation(traces):
    """
    Convert traces to a numerical representation for clustering.
    For simplicity, you can use the length of the trace as a feature.
    """
    return np.array([len(trace) for trace in traces]).reshape(-1, 1)

def find_closest_trace_to_center(traces, centers):
    """
    Find the closest trace to each cluster center.
    """
    closest_traces = []
    for center in centers:
        min_dist = float('inf')
        closest_trace = None
        for trace in traces:
            trace_rep = len(trace)  # Use the length of the trace as the feature
            dist = abs(trace_rep - center)
            if dist < min_dist:
                min_dist = dist
                closest_trace = trace
        closest_traces.append(closest_trace)
    return closest_traces

def perform_kcenter_clustering(y, list_of_vals, k):
    """
    Perform K-center clustering on the Levenshtein distance matrix.
    y: The condensed one-dimensional distance matrix.
    list_of_vals: The list of trace names.
    k: The number of clusters.
    """
    dist_matrix = squareform(y)
    centers = []
    remaining_indices = list(range(len(dist_matrix)))

    # Randomly choose the first center
    first_center = random.choice(remaining_indices)
    centers.append(first_center)
    remaining_indices.remove(first_center)

    while len(centers) < k:
        _, distances_to_closest_center = pairwise_distances_argmin_min(dist_matrix[remaining_indices],
                                                                       dist_matrix[centers])
        next_center = remaining_indices[np.argmax(distances_to_closest_center)]
        centers.append(next_center)
        remaining_indices.remove(next_center)

    return [list_of_vals[center] for center in centers]

def majority_voting(paths):
    """
    Calculate the consensus path using majority voting.

    paths: List of activity sequences, where each sequence is a list of activities.

    Returns:
    consensus_path: The consensus path derived from the majority voting.
    """
    # Determine the maximum length of the paths
    max_length = max(len(path) for path in paths)

    # Pad all paths to the maximum length using "-"
    padded_paths = [path + ['-'] * (max_length - len(path)) for path in paths]

    # Perform majority voting on each position
    consensus_path = []
    for i in range(max_length):
        # Get the activities at the i-th position in all paths
        activities = [path[i] for path in padded_paths]

        # Count occurrences of each activity
        activity_count = Counter(activities)

        # Find the most common activity
        most_common_activity, count = activity_count.most_common(1)[0]

        # If the most common activity appears in more than half of the paths, choose it
        if count > len(paths) // 2:
            consensus_path.append(most_common_activity)
        else:
            # If no activity is a majority, randomly pick one of the non-dash activities
            non_dash_activities = [act for act in activities if act != '-']
            if non_dash_activities:
                consensus_path.append(random.choice(non_dash_activities))
            else:
                consensus_path.append('-')

    return ''.join(consensus_path)


def extract_paths(log):
    """
    Extract the sequences of activities (paths) from the event log.
    """
    paths = []
    for trace in log:
        path = [event['concept:name'] for event in trace]
        paths.append(path)
    return paths

def perform_multi_alignment(filtered_log):
    """
    Perform multi-alignment on the filtered log using the Petri net model.
    """
    # Load Petri net and markings
    pnml_path = "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/example_petri_net2.pnml"
    net, initial_marking, final_marking = petri_importer.apply(pnml_path)

    # Parameters for the multi-alignment
    THETA = 1.1
    MU = 20

    # Perform multi-alignment
    multialignment_result = multii(filtered_log, net, initial_marking, final_marking,
                                   parameters={Parameters.EXPONENT: THETA, Parameters.MARKING_LIMIT: MU})

    # Output the original multi-alignment result
    print("Multi-alignment:", multialignment_result['multi-alignment'])

    # Convert the multi-alignment to a sequence
    sequence = ''.join(
        [activity for transition, activity in multialignment_result['multi-alignment'] if activity is not None])

    # Output the formatted sequence
    print("Aligned Sequence:", sequence)

    print("Maximal Levenshtein Edit Distance to Log:", multialignment_result['max_distance_to_log'])

def print_concept_names(log):
    concept_names = set(trace.attributes["concept:name"] for trace in log)
    print(f"All concept:name values: {concept_names}")


def execute_script():
    log = pm4py.read_xes("/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/Lfull_attribute.xes",
                         return_legacy_log_object=True)
    log = pm4py.sample_cases(log, num_cases=2768)

    print_concept_names(log)

    # Filter traces where concept:name is "Case1", "Case3", or "Case5"
    filtered_log = filter_logs(log, {"Case1", "Case2", "Case3", "Case4", "Case5", "Case7", "Case8", "Case9"})

    # Print the number of traces after filtering
    print(f"Filtered Log Length: {len(filtered_log)}")

    highest_frequency_trace = find_case_with_highest_frequency(filtered_log, {"Case1", "Case2", "Case3", "Case4", "Case5", "Case7", "Case8", "Case9"})
    print(f"Case with the highest frequency of traces: {highest_frequency_trace}")
    # Calculate the Levenshtein distance matrix
    y, list_of_vals = calculate_distance_matrix(filtered_log)

    # Select the new medoid
    if y.size > 0:  # Ensure y is not empty
        new_medoid = find_new_medoid(y, list_of_vals)
        print(f"New Medoid: {new_medoid}")
    else:
        print("Distance matrix is empty, unable to find medoid.")

    # Extract the actual paths (sequences of activities) from the filtered_log
    paths = extract_paths(filtered_log)

    # Calculate the consensus path using majority voting
    consensus_path = majority_voting(paths)
    print(f"Consensus Path: {consensus_path}")

    # Perform K-means to get the closest traces
    trace_representations = get_trace_representation(paths)
    kmeans = KMeans(n_clusters=3, random_state=0).fit(trace_representations)
    centers = kmeans.cluster_centers_.flatten()
    closest_traces = find_closest_trace_to_center(paths, centers)
    print(f"K-means centers: {[''.join(trace) for trace in closest_traces]}")

    # Perform K-center
    k = 3  # Define the number of clusters
    kcenter_centers = perform_kcenter_clustering(y, list_of_vals, k)
    print(f"K-center centers: {kcenter_centers}")

    # Perform multi-alignment
    # log_path = "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/L_3traces.xes"
    pnml_path = "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/ourpetrinet.pnml"
    # log = xes_importer.apply(log_path)
    net, marking, fmarking = petri_importer.apply(pnml_path)

    THETA = 1.1
    MU = 20
    multiali = multii(filtered_log, net, marking, fmarking,
                      parameters={Parameters.EXPONENT: THETA, Parameters.MARKING_LIMIT: MU})
    print("Multi-alignment:", multiali['multi-alignment'])
    print("Maximal Levenshtein Edit Distance to Log:", multiali['max_distance_to_log'])

if __name__ == "__main__":
    execute_script()