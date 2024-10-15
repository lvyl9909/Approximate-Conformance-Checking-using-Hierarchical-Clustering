from xml.etree import ElementTree as ET

# Load and parse the XES file
file_path = '/Users/lyuyilin/Desktop/2024winter/COMP90005/datatsets/OneDrive_1_2024-8-24/BPI2016/BPI2016_Questions.xes'
tree = ET.parse(file_path)
root = tree.getroot()

# Extract all concept:name values for events
concept_names = set()
for event in root.iter('string'):
    if event.attrib.get('key') == 'concept:name':
        concept_names.add(event.attrib.get('value'))

concept_names