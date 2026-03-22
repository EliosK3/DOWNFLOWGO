import os
import csv
import downflowgo.check_files as check_files


class DataManager:
    def __init__(self, config: object):
        self.config = config
        os.makedirs(self.config.path_to_eruptions, exist_ok=True)  # No error raise if keeping existing folder.

        if config.grid_mode == 'yes':
            check_files.overwrite_check_files(config.delete_existing, config.path_to_grid_folder)

        # First check that DEM has a valid format
        check_files.check_dem_valid(config.dem)

    def csv_vent_writer(self, name_vent, easting, northing):
        """
        Create a CSV file containing the name of the vent and its coordinates (latitude and longitude).

        Parameters
        ----------
        name_vent: str
        Name of the vent.

        northing: float
        Latitude.

        easting: float
        Longitude.
        """
        # If csv_vent_file is 0, a new csv file is created based on the given coordinates in the config file
        if self.config.csv_vent_file == '0':
            self.config.from_vent = True
            self.config.csv_vent_file = os.path.join(self.config.path_to_eruptions,
                                                     'csv_vent_file_from_coordinates.csv')
            if os.path.exists(self.config.csv_vent_file):
                os.remove(self.config.csv_vent_file)
            with open(self.config.csv_vent_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter=';')
                writer.writerow(['flow_id', 'X', 'Y'])
                writer.writerow([name_vent, easting, northing])
            print(f"[INFO] New csv created : {self.config.csv_vent_file}")

    def csv_vent_reader(self) -> list:
        """ Read the csv file containing the 'flow_id' and its position 'X' and 'Y' """
        print(f"[INFO] Csv file used : {self.config.csv_vent_file}")
        if not os.path.exists(self.config.csv_vent_file):
            raise FileNotFoundError(f"[ERREUR] File CSV does not exist : {self.config.csv_vent_file}")
        check_files.validate_csv_format(self.config.csv_vent_file)
        with open(self.config.csv_vent_file, 'r') as csvfile:
            csvreader = csv.DictReader(csvfile, delimiter=';')
            return list(csvreader)
