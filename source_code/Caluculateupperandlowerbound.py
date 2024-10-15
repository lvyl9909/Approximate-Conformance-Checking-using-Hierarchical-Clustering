import pm4py
from pm4py.objects.log.importer.xes import importer as xes_importer
from difflib import SequenceMatcher
from collections import Counter

def calculate_edit_distance(trace1, trace2):
    # 使用difflib中的SequenceMatcher计算编辑距离
    return len(trace1) + len(trace2) - 2 * sum(triple.size for triple in SequenceMatcher(None, trace1, trace2).get_matching_blocks())

def calculate_lower_fitness(log1, log2):
    lower_fitness = {}
    trace_count = Counter([trace.attributes["concept:name"] for trace in log1])
    total_weighted_fitness = 0
    total_weight = sum(trace_count.values())  # 总出现次数

    print("Calculating lower fitness for each trace:")
    print("-" * 50)

    for trace_name, count in trace_count.items():
        trace_length = None
        for trace in log1:
            if trace.attributes["concept:name"] == trace_name:
                trace_length = len(trace)
                break

        if trace_length is None:
            continue

        # 计算trace与log2中所有trace的编辑距离并选择最小的
        min_edit_distance = float('inf')
        for model_trace in log2:
            edit_distance = calculate_edit_distance([event['concept:name'] for event in trace], model_trace)
            min_edit_distance = min(min_edit_distance, edit_distance)

        fitness = 1 - min_edit_distance / (trace_length + min(len(x) for x in log2))
        lower_fitness[trace_name] = fitness

        # # print calculation process
        # print(f"Trace Name: {trace_name}")
        # print(f"Edit Distance: {min_edit_distance}")
        # print(f"Trace Length: {trace_length}")
        # print(f"Minimum Length in log2: {min(len(x) for x in log2)}")
        # print(f"Calculated Fitness: {fitness}")
        # print("-" * 50)

        # weighted accumulation
        total_weighted_fitness += fitness * count

    overall_lower_fitness = total_weighted_fitness / total_weight if total_weight > 0 else 0

    return lower_fitness, overall_lower_fitness

def calculate_upper_fitness(log1, log2):

    shortest_path_length = min([len(trace) for trace in log2])

    longest_path_length = max([len(trace) for trace in log2])

    upper_fitness = {}
    trace_count = Counter([trace.attributes["concept:name"] for trace in log1])

    total_weighted_fitness = 0
    total_weight = sum(trace_count.values())  # 总出现次数

    print("Calculating upper fitness for each trace:")
    print("-" * 50)

    for trace_name, count in trace_count.items():
        trace_length = None
        for trace in log1:
            if trace.attributes["concept:name"] == trace_name:
                trace_length = len(trace)
                break

        if trace_length is None:
            continue

        if trace_length < shortest_path_length:
            fitness = 1 - (shortest_path_length - trace_length) / (trace_length + min(len(x) for x in log2))
        elif trace_length > longest_path_length:
            fitness = 1 - (trace_length - longest_path_length) / (trace_length + min(len(x) for x in log2))
        else:
            fitness = 1

        upper_fitness[trace_name] = fitness

        total_weighted_fitness += fitness * count

    overall_upper_fitness = total_weighted_fitness / total_weight if total_weight > 0 else 0

    return upper_fitness, overall_upper_fitness

def main():
    # 导入XES格式的日志文件 log1
    log1_path = "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/L_3traces-added.xes"
    log1 = xes_importer.apply(log1_path)

    # 手动输入log2 trace (log2是多个trace的集合)
    trace1 = ["a", "b", "c", "d", "e", "h"]
    trace2 = ["a", "b", "h"]
    log2 = [trace1, trace2]  # log2 由多个trace组成

    # 计算Fitness下界值
    lower_fitness, overall_lower_fitness = calculate_lower_fitness(log1, log2)
    print("\nLower Fitness for each trace:", lower_fitness)
    print("Overall Lower Fitness:", overall_lower_fitness)
    upper_fitness, overall_upper_fitness = calculate_upper_fitness(log1, log2)
    print("\nUpper Fitness for each trace:", upper_fitness)
    print("Overall Upper Fitness:", overall_upper_fitness)
    final_fitness = (overall_lower_fitness + overall_upper_fitness) / 2
    print("\nFinal Approximate Fitness:", final_fitness)
if __name__ == "__main__":
    main()
