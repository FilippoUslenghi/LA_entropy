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

    def get_LA_map(self) -> ET.Element | None:
        """
        Get the LA map from the study XML file.
        """

        # Check if the file exists
        if not os.path.exists(self.study_path):
            print("File not found.")
            return
        tree = ET.parse(self.study_path)
        # Get the root of the XML file
        root = tree.getroot()

        # Find the 'Maps' tag
        carto_maps = root.find("Maps")

        # Check if the 'Maps' tag exists
        if carto_maps is None:
            print("Tag 'Maps' not found in the XML file.")
            return

        # Check if there are any maps in the study
        if carto_maps.get("Count") == "0":
            print("No maps found.")
            return

        # Get the LA map in the study
        LA_maps_count = 0
        for carto_map in carto_maps.findall("Map"):
            # Get the map name
            map_name = carto_map.get("Name")
            if map_name is None:
                print("Map name not found.")
                continue
            if "LA" in map_name and "Re" not in map_name:
                LA_map = carto_map
                LA_maps_count += 1
        # Check that there is only one LA map
        if LA_maps_count == 0 or LA_maps_count > 1:
            print("There should be only one LA map.")
            return

        # Check if the LA map is not found
        if LA_map is None:
            print("LA map not found.")
            return

        return LA_map

    def get_points(self, map_: ET.Element) -> None:
        """
        Get the points from the LA map.
        """

        # Get the number of points in the LA map
        CartoPoints = map_.find("CartoPoints")
        if CartoPoints is None:
            print("Tag 'CartoPoints' not found in the XML file.")
            return
        points_count = CartoPoints.get("Count")
        if points_count is None:
            print("Number of points not found.")
            return
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
