import csv
import os
import re
import math
import xml.etree.ElementTree as ET
import sys


class Carser:
    """
    A parser for Carto3 V6 data files.
    """

    def __init__(self, path_to_study: str):
        """
        Initialize the Carser object.
        """

        self.study_path = path_to_study
        self.patient = study_path.split("/")[4]
        self.LA_map = None
        self.LA_map_name = None
        self.LA_points = None
        self.LA_mesh_file = None

        self.get_study()

    def get_study(self) -> None:
        """
        Get the study XML file.
        """

        # Check if the file exists
        if not os.path.exists(self.study_path):
            raise FileNotFoundError(f"File {self.study_path} not found.")
        xml_tree = ET.parse(self.study_path)

        try:
            self.LA_map = self.get_LA_map(xml_tree)
        except Exception as e:
            if (
                len(e.args) > 0
                and e.args[0]
                == "The ratio of the number of points between the first and second LA maps is too low."
            ):
                print(f"Skipping patient {self.patient} due to low point ratio.")
                return
            else:
                # Print the state of the object
                print(self.__dict__)
                # Stop the execution of the program
                raise

        self.LA_map_name = self.LA_map.get("Name")
        self.LA_points = self.get_points(self.LA_map)

        self.LA_mesh_file = self.LA_map.get("FileNames")

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
            raise ValueError(
                "No LA map found in the study."
            )  # TODO: some patients have no LA map, skip them?

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

    def get_points(self, map_: ET.Element) -> None:
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
        print(f"Number of points: {points_count}")

        # Get the paths to the Point_Export files
        points_Export_file = os.path.join(
            "/",
            *self.study_path.split(sep="/")[:-1],
            f"{map_.get('Name')}_Points_Export.xml",
        )
        points_Export = ET.parse(points_Export_file)


if __name__ == "__main__":
    data_dir = "/workspace/raw_data/CARTO_DAVIDE/"

    for patient in os.listdir(data_dir):
        if not os.path.isdir(os.path.join(data_dir, patient)):
            continue
        if patient == "Export_REDO-PVI-07_09_2024-13-39-51" or patient == "to-process":
            continue

        subfolder = os.listdir(os.path.join(data_dir, patient))[0]
        study_descriptions = subfolder.replace("Export_", "").split("-")
        study_xml = [
            f
            for f in os.listdir(os.path.join(data_dir, patient, subfolder))
            if all(desc in f for desc in study_descriptions)
        ][0]

        # subfolder = os.listdir(os.path.join(data_dir, patient))[0]
        # study_desc = re.findall(r"[A-Z]+", subfolder.replace("Export_", ""))
        # # Find the study xml file
        # study_xml = [
        #     f
        #     for f in os.listdir(os.path.join(data_dir, patient, subfolder))
        #     if " " in f
        # ]
        # study_xml = [f for f in study_xml if all(desc in f for desc in study_desc) in f]
        # # Create the path to the study
        # study_path = os.path.join(data_dir, patient, subfolder, study_xml)
        study_path = os.path.join(data_dir, patient, subfolder, study_xml)
        carser = Carser(
            path_to_study=study_path,
        )
