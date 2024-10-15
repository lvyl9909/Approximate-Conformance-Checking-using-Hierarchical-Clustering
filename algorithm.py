from scipy.cluster.hierarchy import to_tree, linkage
from pm4py.statistics.attributes.log import get as attributes_filter
import merge_log
import evaluation
from pm4py.objects.conversion.log import converter as log_converter
from enum import Enum
from pm4py.util import exec_utils
from typing import Optional, Dict, Any, Union
from pm4py.objects.log.obj import EventLog, EventStream
from evaluation import eval_raw_leven
import pandas as pd


class Variants(Enum):
    VARIANT_DMM_LEVEN = evaluation.eval_DMM_leven
    VARIANT_AVG_LEVEN = evaluation.eval_avg_leven
    VARIANT_DMM_VEC = evaluation.eval_DMM_variant
    VARIANT_AVG_VEC = evaluation.eval_avg_variant
    DFG = evaluation.dfg_dist
    VARIANT_RAW_LEVEN = evaluation.eval_raw_leven


VARIANT_DMM_LEVEN = Variants.VARIANT_DMM_LEVEN
VARIANT_AVG_LEVEN = Variants.VARIANT_AVG_LEVEN
VARIANT_DMM_VEC = Variants.VARIANT_DMM_VEC
VARIANT_AVG_VEC = Variants.VARIANT_AVG_VEC
DFG = Variants.DFG
VARIANT_RAW_LEVEN = Variants.VARIANT_RAW_LEVEN

VERSIONS = {VARIANT_DMM_LEVEN, VARIANT_AVG_VEC, VARIANT_DMM_VEC, VARIANT_AVG_VEC, DFG, VARIANT_RAW_LEVEN}


def bfs(tree):
    queue = []
    output = []
    queue.append(tree)
    while queue:
        root = queue.pop(0)
        if len(root['children']) > 0:
            name = [root['name']]
            for child in root['children']:
                queue.append(child)
                name.append(child['name'])
            output.append(name)

    return output


def apply(log: Union[EventLog, EventStream, pd.DataFrame], trace_attribute: str, variant=VARIANT_DMM_LEVEN,
          parameters: Optional[Dict[Any, Any]] = None) -> Any:
    if parameters is None:
        parameters = {}

    log = log_converter.apply(log, variant=log_converter.Variants.TO_EVENT_LOG, parameters=parameters)

    percent = 1
    alpha = 0

    list_of_vals = []
    list_log = []
    list_of_vals_dict = attributes_filter.get_trace_attribute_values(log, trace_attribute)

    list_of_vals_keys = list(list_of_vals_dict.keys())
    for i in range(len(list_of_vals_keys)):
        list_of_vals.append(list_of_vals_keys[i])

    for i in range(len(list_of_vals)):
        logsample = merge_log.log2sublog(log, list_of_vals[i], trace_attribute)
        list_log.append(logsample)

    y = exec_utils.get_variant(variant)(list_log, percent, alpha)

    Z = linkage(y, method='average')

    # Create dictionary for labeling nodes by their IDs
    id2name = dict(zip(range(len(list_of_vals)), list_of_vals))

    T = to_tree(Z, rd=False)

    def add_node_with_distance(node, d3Dendro):
        current_node = dict(children=[], name=str(node.id), distance=node.dist, node_id=node.id)
        d3Dendro['children'].append(current_node)

        if node.left:
            add_node_with_distance(node.left, current_node)
        if node.right:
            add_node_with_distance(node.right, current_node)

    d3Dendro = dict(children=[], name="Root1", distance=T.dist)
    add_node_with_distance(T, d3Dendro)

    leafname = merge_log.label_tree(d3Dendro["children"][0], id2name)
    d3Dendro = d3Dendro["children"][0]
    d3Dendro["name"] = 'root'
    tree = d3Dendro

    trilist = bfs(tree)
    trilist[0][0] = trilist[0][1] + '-' + trilist[0][2]

    rootlist = []
    for ele in trilist:
        rootlist.append(ele[0])

    return tree, leafname, Z
