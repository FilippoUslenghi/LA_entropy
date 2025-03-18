import os
import sys
import xml.etree.ElementTree as ET
import csv
import scipy.io as sio
import logging


# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format="%(name)s - %(levelname)s - %(message)s",
    filename="carser.log",
    encoding="utf-8",
    filemode="w",
)
logger = logging.getLogger(__name__)


class Carser:
    """
    A parser for Carto3 V6 data files.
    """

    def __init__(self, path_to_study_xml: str):
        """
        Initialize the Carser object.
        """
        # Check if the file exists
        if not os.path.exists(path_to_study_xml):
            raise FileNotFoundError(f"File {path_to_study_xml} not found.")
        self.study_xml_path = path_to_study_xml
        self.study_dir = os.path.dirname(self.study_xml_path)
        self.patient = study_path.split("/")[4]
        self.LA_map = None
        self.LA_map_name = None
        self.LA_mesh_file = None
        self.LA_mesh = None
        self.LA_points = []
        self.points_signals = {}

    def parse_study(self) -> None:
        """
        Get the study XML file.
        """
        xml_tree = ET.parse(self.study_xml_path)

        try:
            # Get the LA map from the study XML file
            self.LA_map = self.get_LA_map(xml_tree)
        except Exception as e:
            if (
                e.args[0]
                == "The ratio of the number of points between the first and second LA maps is too low."
            ):
                print(f"Skipping patient {self.patient} due to low point ratio.")
                return
            elif e.args[0] == "No LA map found in the study.":
                print(f"Skipping patient {self.patient} due to no LA map.")
                return
            else:
                # Print the state of the object for debugging
                print(self.__dict__)
                # Stop the execution of the program
                raise

        # Get the name of the LA map
        self.LA_map_name = self.LA_map.get("Name")
        # Get the mesh file of the LA map
        self.LA_mesh_file = self.LA_map.get("FileNames")
        if self.LA_mesh_file is None:
            raise ValueError("Attribute 'FileNames' not found in the XML file")
        # Get the mesh of the LA map
        self.LA_mesh = self.get_mesh(self.LA_mesh_file)
        # Get the points from the LA map
        self.LA_points = self.get_points(self.LA_map)
        # Get the data from the points
        for point in self.LA_points:
            self.points_signals[point.get("ID")] = self.get_signals(point)
            return  # Debugging

    __call__ = parse_study

    def get_LA_map(self, xml_tree) -> ET.Element:
        """
        Get the LA map from the study XML file.
        """
        # Get the root of the XML file
        root = xml_tree.getroot()

        # Find the 'Maps' tag
        carto_maps = root.find("Maps")

        # Check if the 'Maps' tag exists
        if carto_maps is None:
            raise ValueError("Tag 'Maps' not found in the XML file.")

        # Check if there are any maps in the study
        if carto_maps.get("Count") == "0":
            raise ValueError("No maps found in the study.")

        # Get the LA map in the study with the highest number of points
        LA_map = None
        LA_maps = []
        bad_strings = ["Re", "POST", "RF"]
        for carto_map in carto_maps.findall("Map"):
            # Get the map name
            map_name = carto_map.get("Name")
            if (
                map_name is not None
                and "LA" in map_name
                and not any(bad_string in map_name for bad_string in bad_strings)
            ):
                n_points = int(carto_map.find("CartoPoints").get("Count"))
                LA_maps.append((carto_map, n_points))

        # Check if there are any LA maps in the study
        if len(LA_maps) == 0:
            raise ValueError("No LA map found in the study.")

        # Sort the LA maps by the number of points
        sorted_LA_maps = sorted(LA_maps, key=lambda x: x[1], reverse=True)
        # Get the first LA map with the highest number of points
        LA_map = sorted_LA_maps[0][0]

        # Check if there are more than one LA map
        if len(sorted_LA_maps) > 1:
            # Check if the ratio of the number of points between the first and second LA maps is too low
            if (
                sorted_LA_maps[0][1] / (sorted_LA_maps[1][1] + sys.float_info.epsilon)
                < 1.5
            ):
                raise ValueError(
                    "The ratio of the number of points between the first and second LA maps is too low."
                )

        # Check if one LA map was found
        if LA_map is None:
            raise ValueError("No LA map found in the study.")

        return LA_map

    def get_points(self, map_: ET.Element) -> list[ET.Element]:
        """
        Get the points from the LA map.
        """
        # Get the number of points in the LA map
        CartoPoints = map_.find("CartoPoints")
        if CartoPoints is None:
            raise ValueError("Tag 'CartoPoints' not found in the XML file")
        points_count = CartoPoints.get("Count")
        if points_count is None:
            raise ValueError("Attribute 'Count' not found in the XML file")

        # Get the points from the LA map
        points_Export_file = os.path.join(
            self.study_dir,
            f"{map_.get('Name')}_Points_Export.xml",
        )
        points_export_collection = ET.parse(points_Export_file).getroot()
        points_list = points_export_collection.findall("Point")

        return points_list

    def get_signals(self, point: ET.Element) -> dict[str, list[str]]:
        # Get the path to the point export file
        point_export_filename = point.get("File_Name")
        if point_export_filename is None:
            raise ValueError("Attribute 'File_Name' not found in the XML file")
        point_export_path = os.path.join(
            self.study_dir,
            point_export_filename,
        )
        point_export_tree = ET.parse(point_export_path)
        point_export_root = point_export_tree.getroot()
        ECG_file = point_export_root.find("ECG")
        if ECG_file is None:
            raise ValueError("Tag 'ECG' not found in the XML file")
        ECG_file_name = ECG_file.get("FileName")
        if ECG_file_name is None:
            raise ValueError("Attribute 'FileName' not found in the XML file")
        ECG_file_path = os.path.join(
            self.study_dir,
            ECG_file_name,
        )

        with open(ECG_file_path, "r") as file:
            # Skip the first two lines
            file.readline()
            file.readline()
            # Read the ECG file header
            header = file.readline().split()
            # Remove the text inside parentheses from the column names
            header = [col.split("(")[0] for col in header]
            # Add an empty string to the header because of formatting issues
            header.append("")
            # Read the ECG file
            reader = csv.reader(
                file,
                delimiter=" ",
                skipinitialspace=True,
                strict=True,
                quoting=csv.QUOTE_NONNUMERIC,
            )
            # Read the data
            data = list(reader)
            # Create a DataFrame
            data_dict = {key: list(val) for key, val in zip(header, zip(*data))}
            # Check if last key is an empty string and remove it
            if data_dict[""]:
                del data_dict[""]

        return data_dict

    def get_electrode_positions(self, point_xml_root: ET.Element):
        pass

    def get_mesh(self, mesh_file: str) -> dict[str, list[tuple]]:
        """
        Get the mesh of the LA map.
        """
        # Get the path to the mesh file
        mesh_file_path = os.path.join(
            self.study_dir,
            mesh_file,
        )
        # Parse the mesh file
        mesh = {}
        vertices = []
        GroupID = []
        triangles = []
        with open(mesh_file_path, "r", encoding="iso-8859-1") as file:
            line = file.readline()
            while "[VerticesSection]" not in line:
                line = file.readline()
            line = file.readline()
            line = file.readline()
            line = file.readline()
            # Load the vertices
            while line != "\n":
                # Split the line by whitespaces
                line_split = line.split()
                # Append the values to the lists
                x = float(line_split[2])
                y = float(line_split[3])
                z = float(line_split[4])
                vertices.append((x, y, z))
                GroupID.append(int(line_split[8]))
                line = file.readline()

            # Load the triangles
            while "[TrianglesSection]" not in line:
                line = file.readline()
            line = file.readline()
            line = file.readline()
            line = file.readline()
            while line != "\n":
                # Split the line by whitespaces
                line_split = line.split()
                # Append the values to the list
                triangles.append(
                    (
                        int(line_split[2]),
                        int(line_split[3]),
                        int(line_split[4]),
                    )
                )
                line = file.readline()

        mesh["vertices"] = vertices
        mesh["triangles"] = triangles
        mesh["GroupID"] = GroupID

        return mesh


if __name__ == "__main__":
    data_dir = "/workspace/raw_data/CARTO_DAVIDE/"

    for patient in os.listdir(data_dir):
        if not os.path.isdir(os.path.join(data_dir, patient)):
            continue
        if patient == "to-process":
            continue
        if patient == "Export_REDO-PVI-07_09_2024-13-39-51":
            carser = Carser(
                os.path.join(data_dir, patient, "REDO PVI 07_09_2024 13-39-51.xml")
            )
            carser()
            continue

        out_dir = os.path.join(data_dir.replace("raw_data", "processed_data"), patient)
        os.makedirs(out_dir, exist_ok=True)
        # Get the subfolders of the patient
        subfolders = os.listdir(os.path.join(data_dir, patient))
        # Check if list is empty
        if not subfolders:
            print(f"Skipping patient {patient} due to no data at all.")
            continue
        subfolder = subfolders[0]
        study_keywords = subfolder.replace("Export_", "").split("-")
        study_xml = [
            f
            for f in os.listdir(os.path.join(data_dir, patient, subfolder))
            if all(desc in f for desc in study_keywords)
        ][0]

        study_path = os.path.join(data_dir, patient, subfolder, study_xml)
        carser = Carser(study_path)
        logging.debug(f"Processing patient {patient}")
        carser()
        sio.savemat(
            os.path.join(out_dir, "LA_mesh.mat"),
            carser.LA_mesh,
            oned_as="column",
        )
        sio.savemat(
            os.path.join(out_dir, "points_signals.mat"),
            carser.points_signals,
            oned_as="column",
        )
        break  # Debugging

        # TODO: use numpy to store the data in a more efficient way
        # TODO: maybe use pandas to store the data in a more efficient way
