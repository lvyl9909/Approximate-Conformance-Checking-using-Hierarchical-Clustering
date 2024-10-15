
from pm4py.objects.petri_net.obj import PetriNet, Marking
from pm4py.objects.petri_net.utils import petri_utils
from pm4py.objects.petri_net.exporter import exporter as pnml_exporter
from pm4py.visualization.petri_net import visualizer as pn_visualizer

# 1. 初始化一个新的Petri网
net = PetriNet("complex_petri_net_with_invisible")

# 2. 创建位置
places = {
    "p_start": PetriNet.Place("p_start"),
    "p1": PetriNet.Place("p1"),
    "p2": PetriNet.Place("p2"),
    "p3": PetriNet.Place("p3"),
    "p4": PetriNet.Place("p4"),
    "p5": PetriNet.Place("p5"),
    "p6": PetriNet.Place("p6"),
    "p7": PetriNet.Place("p7"),
    "p8": PetriNet.Place("p8"),
    "p9": PetriNet.Place("p9"),
    "p10": PetriNet.Place("p10"),
    "p_end": PetriNet.Place("p_end")
}

# 将位置添加到Petri网
for place in places.values():
    net.places.add(place)

# 3. 创建变迁（包括可见和不可见的）
transitions = {
    "t_a": PetriNet.Transition("t_a", "a"),
    "t_b": PetriNet.Transition("t_b", "b"),
    "t_c": PetriNet.Transition("t_c", "c"),
    "t_d": PetriNet.Transition("t_d", "d"),
    "t_e": PetriNet.Transition("t_e", "e"),
    "t_f": PetriNet.Transition("t_f", "f"),
    "t_g": PetriNet.Transition("t_g", "g"),
    "t_h": PetriNet.Transition("t_h", "h"),
    "t_invisible1": PetriNet.Transition("t_invisible1", None),  # invisible transition
    "t_invisible2": PetriNet.Transition("t_invisible2", None),  # invisible transition
    "t_invisible3": PetriNet.Transition("t_invisible3", None),  # invisible transition
    "t_invisible4": PetriNet.Transition("t_invisible4", None),  # invisible transition
    "t_invisible5": PetriNet.Transition("t_invisible5", None),  # invisible transition
    "t_invisible6": PetriNet.Transition("t_invisible6", None),  # invisible transition
    "t_invisible7": PetriNet.Transition("t_invisible7", None)   # invisible transition
}

# 将变迁添加到Petri网
for transition in transitions.values():
    net.transitions.add(transition)

# 4. 连接位置和变迁（包括不可见的变迁）
petri_utils.add_arc_from_to(places["p_start"], transitions["t_a"], net)
petri_utils.add_arc_from_to(transitions["t_a"], places["p1"], net)

petri_utils.add_arc_from_to(places["p1"], transitions["t_invisible1"], net)  # invisible
petri_utils.add_arc_from_to(transitions["t_invisible1"], places["p2"], net)  # invisible

petri_utils.add_arc_from_to(places["p1"], transitions["t_b"], net)  # invisible
petri_utils.add_arc_from_to(transitions["t_b"], places["p2"], net)  # invisible

petri_utils.add_arc_from_to(places["p2"], transitions["t_c"], net)
petri_utils.add_arc_from_to(transitions["t_c"], places["p3"], net)

petri_utils.add_arc_from_to(places["p2"], transitions["t_invisible2"], net)
petri_utils.add_arc_from_to(transitions["t_invisible2"], places["p3"], net)

petri_utils.add_arc_from_to(places["p3"], transitions["t_d"], net)
petri_utils.add_arc_from_to(transitions["t_d"], places["p4"], net)

petri_utils.add_arc_from_to(places["p3"], transitions["t_invisible3"], net)
petri_utils.add_arc_from_to(transitions["t_invisible3"], places["p4"], net)

petri_utils.add_arc_from_to(places["p4"], transitions["t_invisible4"], net)
petri_utils.add_arc_from_to(transitions["t_invisible4"], places["p10"], net)

petri_utils.add_arc_from_to(places["p4"], transitions["t_invisible5"], net)
petri_utils.add_arc_from_to(transitions["t_invisible5"], places["p5"], net)
petri_utils.add_arc_from_to(transitions["t_invisible5"], places["p6"], net)

petri_utils.add_arc_from_to(places["p5"], transitions["t_f"], net)
petri_utils.add_arc_from_to(transitions["t_f"], places["p7"], net)

petri_utils.add_arc_from_to(places["p7"], transitions["t_invisible6"], net)
petri_utils.add_arc_from_to(transitions["t_invisible6"], places["p10"], net)

petri_utils.add_arc_from_to(places["p6"], transitions["t_e"], net)
petri_utils.add_arc_from_to(transitions["t_e"], places["p8"], net)

petri_utils.add_arc_from_to(places["p6"], transitions["t_invisible7"], net)
petri_utils.add_arc_from_to(transitions["t_invisible7"], places["p9"], net)

petri_utils.add_arc_from_to(places["p8"], transitions["t_g"], net)
petri_utils.add_arc_from_to(transitions["t_g"], places["p9"], net)

petri_utils.add_arc_from_to(places["p9"], transitions["t_invisible6"], net)


petri_utils.add_arc_from_to(places["p10"], transitions["t_h"], net)
petri_utils.add_arc_from_to(transitions["t_h"], places["p_end"], net)
# 5. 定义初始标识和最终标识
initial_marking = Marking()
initial_marking[places["p_start"]] = 1
final_marking = Marking()
final_marking[places["p_end"]] = 1

# 6. 导出Petri网为PNML文件
pnml_exporter.apply(net, initial_marking, "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/example_petri_net2.pnml", final_marking=final_marking)

# 7. 可选：可视化Petri网
parameters = {pn_visualizer.Variants.WO_DECORATION.value.Parameters.FORMAT: "svg"}
gviz = pn_visualizer.apply(net, initial_marking, final_marking, parameters=parameters)
pn_visualizer.view(gviz)



from pm4py.objects.petri_net.obj import PetriNet, Marking
from pm4py.objects.petri_net.utils import petri_utils
from pm4py.objects.petri_net.exporter import exporter as pnml_exporter
from pm4py.visualization.petri_net import visualizer as pn_visualizer

# 1. 初始化一个新的Petri网
net = PetriNet("complex_petri_net_with_invisible")

# 2. 创建位置
places = {
    "p_start": PetriNet.Place("p_start"),
    "p1": PetriNet.Place("p1"),
    "p2": PetriNet.Place("p2"),
    "p3": PetriNet.Place("p3"),
    "p4": PetriNet.Place("p4"),
    "p5": PetriNet.Place("p5"),
    "p6": PetriNet.Place("p6"),
    "p7": PetriNet.Place("p7"),
    "p8": PetriNet.Place("p8"),
    "p9": PetriNet.Place("p9"),
    "p10": PetriNet.Place("p10"),
    "p_end": PetriNet.Place("p_end")
}

# 将位置添加到Petri网
for place in places.values():
    net.places.add(place)

# 3. 创建变迁（包括可见和不可见的）
transitions = {
    "t_a": PetriNet.Transition("t_a", "a"),
    "t_b": PetriNet.Transition("t_b", "b"),
    "t_c": PetriNet.Transition("t_c", "c"),
    "t_d": PetriNet.Transition("t_d", "d"),
    "t_e": PetriNet.Transition("t_e", "e"),
    "t_f": PetriNet.Transition("t_f", "f"),
    "t_g": PetriNet.Transition("t_g", "g"),
    "t_h": PetriNet.Transition("t_h", "h"),
    "t_invisible1": PetriNet.Transition("t_invisible1", None),  # invisible transition
    "t_invisible2": PetriNet.Transition("t_invisible2", None),  # invisible transition
    "t_invisible3": PetriNet.Transition("t_invisible3", None),  # invisible transition
    "t_invisible4": PetriNet.Transition("t_invisible4", None),  # invisible transition
    "t_invisible5": PetriNet.Transition("t_invisible5", None),  # invisible transition
    "t_invisible6": PetriNet.Transition("t_invisible6", None),  # invisible transition
    # "t_invisible7": PetriNet.Transition("t_invisible7", None)   # invisible transition
}

# 将变迁添加到Petri网
for transition in transitions.values():
    net.transitions.add(transition)

# 4. 连接位置和变迁（包括不可见的变迁）
petri_utils.add_arc_from_to(places["p_start"], transitions["t_a"], net)
petri_utils.add_arc_from_to(transitions["t_a"], places["p1"], net)

# petri_utils.add_arc_from_to(places["p1"], transitions["t_invisible1"], net)  # invisible
# petri_utils.add_arc_from_to(transitions["t_invisible1"], places["p2"], net)  # invisible

petri_utils.add_arc_from_to(places["p1"], transitions["t_b"], net)  # invisible
petri_utils.add_arc_from_to(transitions["t_b"], places["p2"], net)  # invisible

petri_utils.add_arc_from_to(places["p2"], transitions["t_c"], net)
petri_utils.add_arc_from_to(transitions["t_c"], places["p3"], net)

petri_utils.add_arc_from_to(places["p2"], transitions["t_invisible1"], net)
petri_utils.add_arc_from_to(transitions["t_invisible1"], places["p3"], net)

petri_utils.add_arc_from_to(places["p3"], transitions["t_d"], net)
petri_utils.add_arc_from_to(transitions["t_d"], places["p4"], net)

petri_utils.add_arc_from_to(places["p3"], transitions["t_invisible2"], net)
petri_utils.add_arc_from_to(transitions["t_invisible2"], places["p4"], net)

petri_utils.add_arc_from_to(places["p4"], transitions["t_invisible3"], net)
petri_utils.add_arc_from_to(transitions["t_invisible3"], places["p10"], net)

petri_utils.add_arc_from_to(places["p4"], transitions["t_invisible4"], net)
petri_utils.add_arc_from_to(transitions["t_invisible4"], places["p5"], net)
petri_utils.add_arc_from_to(transitions["t_invisible4"], places["p6"], net)

petri_utils.add_arc_from_to(places["p5"], transitions["t_f"], net)
petri_utils.add_arc_from_to(transitions["t_f"], places["p7"], net)

petri_utils.add_arc_from_to(places["p7"], transitions["t_invisible5"], net)
petri_utils.add_arc_from_to(transitions["t_invisible5"], places["p10"], net)

petri_utils.add_arc_from_to(places["p6"], transitions["t_e"], net)
petri_utils.add_arc_from_to(transitions["t_e"], places["p8"], net)

petri_utils.add_arc_from_to(places["p6"], transitions["t_invisible6"], net)
petri_utils.add_arc_from_to(transitions["t_invisible6"], places["p9"], net)

petri_utils.add_arc_from_to(places["p8"], transitions["t_g"], net)
petri_utils.add_arc_from_to(transitions["t_g"], places["p9"], net)

petri_utils.add_arc_from_to(places["p9"], transitions["t_invisible5"], net)


petri_utils.add_arc_from_to(places["p10"], transitions["t_h"], net)
petri_utils.add_arc_from_to(transitions["t_h"], places["p_end"], net)
# 5. 定义初始标识和最终标识
initial_marking = Marking()
initial_marking[places["p_start"]] = 1
final_marking = Marking()
final_marking[places["p_end"]] = 1

# 6. 导出Petri网为PNML文件
pnml_exporter.apply(net, initial_marking, "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/example_petri_net2.pnml", final_marking=final_marking)

# 7. 可选：可视化Petri网
parameters = {pn_visualizer.Variants.WO_DECORATION.value.Parameters.FORMAT: "svg"}
gviz = pn_visualizer.apply(net, initial_marking, final_marking, parameters=parameters)
pn_visualizer.view(gviz)