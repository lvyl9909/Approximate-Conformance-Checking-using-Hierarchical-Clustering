
import xml.etree.ElementTree as ET

def clean_traces(xml_input_path, xml_output_path):
    # 解析 XML 文件
    tree = ET.parse(xml_input_path)
    root = tree.getroot()

    # 遍历每一个 trace
    for trace in root.findall('trace'):
        # 获取 trace 中的所有 string 和 date 标签
        for child in list(trace):
            # 如果该标签是 string 并且 key 属性不是 "concept:name"，则移除该标签
            if child.tag == 'string' and child.attrib.get('key') != 'concept:name':
                trace.remove(child)
            # 如果该标签是 date 或 float，也进行移除
            elif child.tag in ['date', 'float']:
                trace.remove(child)

    # 保存修改后的 XML 文件
    tree.write(xml_output_path, encoding='utf-8', xml_declaration=True)
    print(f"文件处理完成，保存至: {xml_output_path}")

# 输入和输出文件路径
xml_input_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2015/BPIC15_1.xes'
xml_output_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2015/BPIC15.xes'

# 调用函数
clean_traces(xml_input_path, xml_output_path)
