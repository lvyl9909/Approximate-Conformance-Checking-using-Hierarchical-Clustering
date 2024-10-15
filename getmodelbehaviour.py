
#######Get the model behaviour########
import pm4py
from pm4py.objects.log.obj import Trace, Event
from pm4py.objects.petri_net.importer import importer as pnml_importer
from pm4py.algo.conformance.alignments.petri_net import algorithm as alignments

trace = Trace()
trace.append(Event({"concept:name": "a"}))
trace.append(Event({"concept:name": "b"}))
trace.append(Event({"concept:name": "c"}))
trace.append(Event({"concept:name": "d"}))
# trace.append(Event({"concept:name": "e"}))
trace.append(Event({"concept:name": "f"}))
trace.append(Event({"concept:name": "e"}))

# trace.append(Event({"concept:name": "f"}))


pnml_path = "/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/wllpetrinet.pnml"
net, initial_marking, final_marking = pnml_importer.apply(pnml_path)

alignment = alignments.apply(trace, net, initial_marking, final_marking)

print("Optimal-alignment trace:")
for step in alignment["alignment"]:
    print(step)
