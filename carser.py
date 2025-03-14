import csv
import os
import xml.etree.ElementTree as ET


class Carser:
    """
    A parser for Carto3 V6 data files.
    """

    def __init__(self, path_to_study: str):
        """
        Initialize the Carser object.
        """

        self.study_path = path_to_study
        self.LA_map = None
        self.LA_map_name = None
        self.LA_points = None
        self.LA_mesh_file = None

    def get_study(self) -> None:
        """
        Get the study XML file.
        """

        self.LA_map = self.get_LA_map()
        if self.LA_map is None:
            return
        self.LA_map_name = self.LA_map.get("Name")
        self.LA_points = self.get_points(self.LA_map)

        self.LA_mesh_file = self.LA_map.get("FileNames")

    def get_LA_map(self) -> ET.Element:
        """
        Get the LA map from the study XML file.
        """

        # Check if the file exists
        if not os.path.exists(self.study_path):
            raise FileNotFoundError(f"File {self.study_path} not found.")
        tree = ET.parse(self.study_path)
        # Get the root of the XML file
        root = tree.getroot()

        # Find the 'Maps' tag
        carto_maps = root.find("Maps")

        # Check if the 'Maps' tag exists
        if carto_maps is None:
            raise ValueError("Tag 'Maps' not found in the XML file.")

        # Check if there are any maps in the study
        if carto_maps.get("Count") == "0":
            raise ValueError("No maps found in the study.")

        # Get the LA map in the study
        LA_map = None
        LA_maps_count = 0
        for carto_map in carto_maps.findall("Map"):
            # Get the map name
            map_name = carto_map.get("Name")
            if map_name is not None and "LA" in map_name and "Re" not in map_name:
                LA_map = carto_map
                LA_maps_count += 1

        # Check if one LA map was found
        if LA_map is None:
            raise ValueError("No LA map found in the study.")

        # Check that there is only one LA map
        if LA_maps_count == 0 or LA_maps_count > 1:
            raise ValueError("There should be only one LA map in the study.")

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
        print(points_Export)


if __name__ == "__main__":
    carser = Carser(
        "/workspace/raw_data/CARTO_DAVIDE/100/Export_FAREDO-11_03_2023-14-15-05/FAREDO 11_03_2023 14-15-05.xml"
    )
    carser.get_study()
