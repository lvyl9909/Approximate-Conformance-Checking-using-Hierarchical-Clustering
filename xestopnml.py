import pm4py
from pm4py.objects.log.importer.xes import importer as xes_importer
from pm4py.algo.discovery.inductive import algorithm as inductive_miner
from pm4py.objects.petri_net.exporter import exporter as pnml_exporter
from pm4py.objects.conversion.process_tree import converter as pt_converter

# 1. 导入事件日志
log = xes_importer.apply('/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2016/BPI2016_Questions.xes')

# 2. 使用 Inductive Miner - Infrequent 算法生成 Process Tree，设置噪声参数为 0.4
parameters = {"noiseThreshold": 0.4}  # 使用字符串键设置噪声阈值
process_tree = inductive_miner.apply(log, variant=inductive_miner.Variants.IMf, parameters=parameters)

# 3. 将 Process Tree 转换为 Petri 网
net, initial_marking, final_marking = pt_converter.apply(process_tree)

# 4. 将 Petri 网导出为 PNML 格式
pnml_exporter.apply(net, initial_marking, '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2016/BPI2016_Questions.pnml', final_marking=final_marking)

print("Petri net has been exported to PNML format with noise threshold set to 0.4.")
