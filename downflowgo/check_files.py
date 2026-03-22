import csv
import numpy as np
import os
import shutil

def check_dem_valid(dem:str):
    """ Validate that a DEM (ASCII Grid) file is well-formed and internally consistent. """
    expected_headers = ['ncols', 'nrows', 'xllcorner', 'yllcorner', 'cellsize', 'nodata_value']
    with open(dem, 'r') as file:
        # Read and check header
        header_lines = [next(file).strip() for _ in range(6)]
        header_values = {}

        for i, line in enumerate(header_lines):
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"Malformed header line {i + 1} in DEM '{dem}': {line}")

            key, value = parts[0].lower(), parts[1]
            if key != expected_headers[i]:
                raise ValueError(
                    f"Unexpected header key at line {i + 1} in DEM '{dem}': got '{key}', expected '{expected_headers[i]}'"
                )

            header_values[key] = value

        # Convert header values to appropriate types
        try:
            ncols = int(header_values['ncols'])
            nrows = int(header_values['nrows'])
            xllcorner = float(header_values['xllcorner'])
            yllcorner = float(header_values['yllcorner'])
            cellsize = float(header_values['cellsize'])
            nodata_value = float(header_values['nodata_value'])
        except ValueError as e:
            raise ValueError(f"Invalid header value in DEM '{dem}': {e}")


        # Read the remaining lines (data grid)
        data_lines = [line.strip().split() for line in file if line.strip()]

        if len(data_lines) != nrows:
            raise ValueError(f"Expected {nrows} rows of data, but found {len(data_lines)} in DEM '{dem}'.")

        for i, row in enumerate(data_lines):
            if len(row) != ncols:
                raise ValueError(f"Expected {ncols} columns in row {i + 1}, but got {len(row)} in DEM '{dem}'.")

        try:
            data = np.array(data_lines, dtype=float)
        except ValueError:
            raise ValueError(f"Non-numeric value found in data grid of DEM '{dem}'.")

    print("DEM is valid.")
    return True

def check_vent_in_dem(long:float, lat:float, dem:str):
    """ to check dem headers and data lines as well as vent position within the dem """

    long = float(long)
    lat = float(lat)

    with open(dem) as file:
        header_lines = [next(file) for _ in range(6)]
        ncols = int(header_lines[0].split()[1])
        nrows = int(header_lines[1].split()[1])
        xllcorner = float(header_lines[2].split()[1])
        yllcorner = float(header_lines[3].split()[1])
        cellsize = float(header_lines[4].split()[1])

        # Check that (long, lat) is inside the DEM extent
        x_max = xllcorner + ncols * cellsize
        y_max = yllcorner + nrows * cellsize

        if not (xllcorner <= long <= x_max and yllcorner <= lat <= y_max):
            raise ValueError(
                f"The coordinates of the vent (long={long}, lat={lat}) is outside the DEM extent:\n"
                f"x: [{xllcorner}, {x_max}], y: [{yllcorner}, {y_max}]")
    return True

def validate_csv_format(csv_vent_file:str) -> bool:
    expected_header = ['flow_id', 'X', 'Y']

    try:
        with open(csv_vent_file, newline='', encoding='utf-8') as csvfile:
            reader = csv.reader(csvfile, delimiter=';')
            rows = list(reader)

            # Check if file is empty
            if not rows:
                print("The file is empty.")
                return False

            # Check header
            header = rows[0]
            if header != expected_header:
                print(f"Incorrect header. Expected: {expected_header}, Found: {header}")
                return False

            # Check each data row
            for i, row in enumerate(rows[1:], start=2):
                if len(row) != 3:
                    print(f"Row {i} does not have exactly 3 fields: {row}")
                    return False
                flow_id, x_str, y_str = row
                try:
                    float(x_str)
                    float(y_str)
                except ValueError:
                    print(f"Row {i} has invalid float values: {x_str}, {y_str}")
                    return False

        print("CSV format is valid.")
        return True

    except FileNotFoundError:
        print(f"File not found: {csv_vent_file}")
        return False
    except Exception as e:
        print(f"Error while reading file: {e}")
        return False

def overwrite_check_files(delete_existing:bool, path_to_folder:str):
    """ Check if the target folder 'path_to_folder' should be deleted, and create an empty one if it does not exist. """
    if os.path.exists(path_to_folder):
        if delete_existing:
            shutil.rmtree(path_to_folder)
        else:
            answer = input(f"The folder '{path_to_folder}' already exists. Overwrite it? [y/N]: ").strip().lower()
            if answer == "y":
                shutil.rmtree(path_to_folder)
            else:
                print("Keeping existing folder.")

    # Create folder (fresh or because it didn’t exist before)
    os.makedirs(path_to_folder, exist_ok=True) # No error raise if keeping existing folder.