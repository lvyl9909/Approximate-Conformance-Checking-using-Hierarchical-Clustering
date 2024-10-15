
import pm4py
import time
from pm4py.objects.log.importer.xes import importer as xes_importer
from pm4py.objects.petri_net.importer import importer as pnml_importer
from pm4py.algo.conformance.alignments.petri_net import algorithm as alignments
from pm4py.algo.evaluation.replay_fitness import algorithm as replay_fitness

log_path = "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis/Sepsis.xes"
log = xes_importer.apply(log_path)

pnml_path = "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/Sepsis/spesis_0.4.pnml"
net, initial_marking, final_marking = pnml_importer.apply(pnml_path)

start_time_alignments = time.time()

alignments_result = alignments.apply(log, net, initial_marking, final_marking)

end_time_alignments = time.time()

alignment_execution_time = end_time_alignments - start_time_alignments

start_time_fitness = time.time()

fitness = replay_fitness.apply(log, net, initial_marking, final_marking,
                               variant=replay_fitness.Variants.ALIGNMENT_BASED)

end_time_fitness = time.time()

fitness_execution_time = end_time_fitness - start_time_fitness

print("Alignments Result:")
for trace_alignment in alignments_result:
    print(trace_alignment)

print("Fitness:")
print(fitness)

print(f"Alignment Execution Time: {alignment_execution_time:.4f} seconds")
print(f"Fitness Calculation Time: {fitness_execution_time:.4f} seconds")
