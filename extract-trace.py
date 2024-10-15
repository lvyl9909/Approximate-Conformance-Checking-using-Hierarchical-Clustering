import xml.etree.ElementTree as ET


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


if __name__ == "__main__":
    # 解析txt文件中的轨迹代码
    txt_file_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPIC2013/cluster-10.txt'
    center_tracks, most_frequent_tracks = parse_tracks_from_txt(txt_file_path)

    # 从XML文件中提取相应的trace节点并生成新的XML文件
    xml_input_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPIC2013/BPIC2013_inc.xes'
    xml_output_center_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPIC2013/incluster-medoid.xes'
    xml_output_frequent_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPIC2013/incluster-frequency.xes'

    extract_traces_to_xml(
        xml_input_path,
        xml_output_center_path,
        xml_output_frequent_path,
        center_tracks,
        most_frequent_tracks
    )