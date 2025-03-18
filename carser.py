import os
import sys
import xml.etree.ElementTree as ET
import csv
import scipy.io as sio
import logging
import numpy as np
import pandas as pd

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
        self.points_signals["signals"] = []
        self.points_signals["points_IDs"] = []
        for point in self.LA_points:
            logger.debug(
                f"Processing point {point.get('ID')} of patient {self.patient}"
            )
            signals, columns = self.get_signals(point)
            if signals is None:
                continue
            self.points_signals["signals"].append(signals)
            self.points_signals["points_IDs"].append(point.get("ID"))
            self.points_signals["columns"] = columns
            # return  # Debugging

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

    def get_signals(
        self, point: ET.Element
    ) -> tuple[np.ndarray | None, list[str] | None]:
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

        # Check which connectors were used
        connectors = point_export_root.find("Positions").findall("Connector")  # type: ignore
        has_pole_a = False
        has_mcc_dx = False
        if not (has_pole_a or has_mcc_dx):
            for connector in connectors:
                if connector.get("MAGNETIC_20_POLE_A_CONNECTOR") is not None:
                    has_pole_a = True

                if connector.get("MCC_DX_CONNECTOR") is not None:
                    has_mcc_dx = True

            if has_mcc_dx and has_pole_a:
                raise ValueError(
                    "MCC-DX and 20 Pole A connectors present in the same point export."
                )

        if not (has_pole_a or has_mcc_dx):
            print(
                f"Skipping point {point.get('ID')} of patient {self.patient} due to no connectors."
            )
            return (None, None)

        # Get the ECG file
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

        # Get the ECG data
        ECG_columns = [
            "V1(22)",
            "V2(23)",
            "V3(24)",
            "V4(25)",
            "V5(26)",
            "V6(27)",
            "I(110)",
            "II(111)",
            "III(112)",
            "aVL(171)",
            "aVR(172)",
            "aVF(173)",
        ]
        CS_columns = ["CS1-CS2(101)", "CS5-CS6(105)", "CS9-CS10(109)"]
        pole_columns = []
        if has_pole_a:
            pole_columns = [
                "20A_1-2(113)",
                "20A_2-3(114)",
                "20A_3-4(115)",
                "20A_5-6(117)",
                "20A_6-7(118)",
                "20A_7-8(119)",
                "20A_9-10(121)",
                "20A_10-11(122)",
                "20A_11-12(123)",
                "20A_13-14(125)",
                "20A_14-15(126)",
                "20A_15-16(127)",
                "20A_17-18(129)",
                "20A_18-19(130)",
                "20A_19-20(131)",
            ]
        elif has_mcc_dx:
            pole_columns = [
                "MCC_Dx_BiPolar_1(197)",
                "MCC_Dx_BiPolar_2(198)",
                "MCC_Dx_BiPolar_3(199)",
                "MCC_Dx_BiPolar_4(200)",
                "MCC_Dx_BiPolar_5(201)",
                "MCC_Dx_BiPolar_6(202)",
                "MCC_Dx_BiPolar_7(203)",
                "MCC_Dx_BiPolar_8(204)",
                "MCC_Dx_BiPolar_9(205)",
                "MCC_Dx_BiPolar_10(206)",
                "MCC_Dx_BiPolar_11(207)",
                "MCC_Dx_BiPolar_12(208)",
                "MCC_Dx_BiPolar_13(209)",
                "MCC_Dx_BiPolar_14(210)",
                "MCC_Dx_BiPolar_15(211)",
                "MCC_Dx_BiPolar_16(212)",
                "MCC_Dx_BiPolar_17(213)",
                "MCC_Dx_BiPolar_18(214)",
                "MCC_Dx_BiPolar_19(215)",
                "MCC_Dx_BiPolar_20(216)",
                "MCC_Dx_BiPolar_21(217)",
                "MCC_Dx_BiPolar_22(218)",
                "MCC_Dx_BiPolar_23(219)",
                "MCC_Dx_BiPolar_24(220)",
                "MCC_Dx_BiPolar_25(221)",
                "MCC_Dx_BiPolar_26(222)",
                "MCC_Dx_BiPolar_27(223)",
                "MCC_Dx_BiPolar_28(224)",
                "MCC_Dx_BiPolar_29(225)",
                "MCC_Dx_BiPolar_30(226)",
                "MCC_Dx_BiPolar_31(227)",
                "MCC_Dx_BiPolar_32(228)",
                "MCC_Dx_BiPolar_33(229)",
                "MCC_Dx_BiPolar_34(230)",
                "MCC_Dx_BiPolar_35(231)",
                "MCC_Dx_BiPolar_36(232)",
                "MCC_Dx_BiPolar_37(233)",
                "MCC_Dx_BiPolar_38(234)",
                "MCC_Dx_BiPolar_39(235)",
                "MCC_Dx_BiPolar_40(236)",
            ]

        columns = pole_columns + ECG_columns + CS_columns
        signal_data = (
            pd.read_csv(
                ECG_file_path,
                skiprows=2,
                sep=" ",
                skipinitialspace=True,
                usecols=columns,
            )
            .sort_values(axis=0, by=columns)
            .to_numpy(dtype=np.float64)
        )

        # Get the gain value
        with open(ECG_file_path, "r") as file:
            # Skip the first two lines
            file.readline()
            line = file.readline()
            # Get the gain value
            gain = float(line.split("=")[1])

        signal_data = signal_data * gain
        columns = [column.split("(")[0] for column in columns]

        # data_dict = np.empty([1 + 2500, len(columns)], dtype=object)
        # data_dict[0, :] = columns
        # data_dict[1:, :] = signal_data
        # for i, column in enumerate(columns, start=1):
        #     column = column.split("(")[0]
        #     data_dict[column] = signal_data[:, i]

        return signal_data, columns

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
        # break  # Debugging
