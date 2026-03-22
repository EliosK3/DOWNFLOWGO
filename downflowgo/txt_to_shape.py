import fiona
import pandas as pd
import csv
import rasterio
import numpy as np
from rasterio.transform import from_origin
from rasterio.crs import CRS
from shapely.geometry import shape, mapping, LineString, Point

def get_path_shp(losd_file, shp_losd_file, epsg_code):
    # import points from losd_file
    lineDf = pd.read_csv(losd_file, header=0, sep='\t')
    lineDf.head()
    # define schema for line shape file
    schema = {'geometry': 'LineString', 'properties': [('L', 'str')]}
    # open a fiona object and write shp_losd_file.shp
    lineShp = fiona.open(shp_losd_file, mode='w',
                         driver='ESRI Shapefile', schema=schema, crs=f'epsg:{epsg_code}')

    # get list of points
    xyList = []
    rowName = ''
    for index, row in lineDf.iterrows():
        xyList.append((row.x, row.y))
        rowName = row.L

    # save record and close shapefile
    rowDict = {'geometry': {'type': 'LineString', 'coordinates': xyList}, 'properties': {'L': rowName}, }
    lineShp.write(rowDict)
    # close fiona object
    lineShp.close()
    print(f"++++++++++++++++++ Losd file is saved in shape format at: '{shp_losd_file}'+++++++++++++++++")

def get_runouts_shp(run_outs_file, shp_runouts, epsg_code):

    # import points from slope file
    pointDf = pd.read_csv(run_outs_file, header=0, sep=',')
    pointDf.head()
    # define schema for line shape file
    schema = {'geometry': 'Point', 'properties': [("flow_id", 'str'), ("Eff_r", 'int'), ("X_run_out", 'float'),
                                                  ("Y_run_out", 'float'), ("Channel_Depth", 'float'), 
                                                  ("Channel_Width_init", 'float'), ("Elevation_run_out", 'int'), 
                                                  ("Distance_run_out", 'int') ]}
    # open a fiona object
    pointShp = fiona.open(shp_runouts, mode='w',
                         driver='ESRI Shapefile', schema=schema, crs=f'epsg:{epsg_code}')

    # iterate over each row in the dataframe and save record
    for index, row in pointDf.iterrows():
        rowDict = {
            'geometry': {'type': 'Point',
                         'coordinates': (row.X_run_out, row.Y_run_out)},
            'properties': {'flow_id': row.flow_id, 'Eff_r': row.Effusion_rate, 'X_run_out': row.X_run_out,
                           'Y_run_out': row.Y_run_out, 'Channel_Depth': row.Depth, 'Channel_Width_init': row.Width_init,
                           'Elevation_run_out': row.Elevation_run_out, 'Distance_run_out': row.Distance_run_out
                           }}
        pointShp.write(rowDict)
    # close fiona object
    pointShp.close()
    print(f"----------- Runouts coordinates are saved as shape file in:'{shp_runouts}'---------------")

def get_vent_shp(csv_vent_file, shp_vent_file, epsg_code):
    # import points from slope file
    with open(csv_vent_file,'r') as file:
        try:
            dialect = csv.Sniffer().sniff(file.read(1024))
            file.seek(0)  # Reset file pointer after sniffing
        except csv.Error:
            print("Could not determine delimiter.")
            return
    if dialect.delimiter == ';':
        pointDf = pd.read_csv(csv_vent_file, header=0, sep=';')
    elif dialect.delimiter == ',':
        pointDf = pd.read_csv(csv_vent_file, header=0, sep=',')

    pointDf.head()
    # define schema for line shape file
    schema = {'geometry': 'Point', 'properties': [("flow_id", 'str')]}
    # open a fiona object

    pointShp = fiona.open(shp_vent_file,
                          mode='w',
                         driver='ESRI Shapefile',
                          schema=schema,
                          crs=f'epsg:{epsg_code}')

    # iterate over each row in the dataframe and save record
    first_coordinates = None
    for index, row in pointDf.iterrows():
        # Check if 'X_init' and 'Y_init' exist in the row
        if 'X_init' in row and 'Y_init' in row:
            coordinates = (row.X_init, row.Y_init)
        else:
            coordinates = (row.X, row.Y)

        # If it's the first iteration, save the coordinates
        if first_coordinates is None:
            first_coordinates = coordinates

        # Skip the row if it has the same coordinates as the first one
        if coordinates == first_coordinates and index > 0:
            continue  # Skip writing the duplicate points

        rowDict = {
            'geometry': {'type': 'Point',
                         'coordinates': coordinates},
            'properties': {'flow_id': row.flow_id
                           }}
        pointShp.write(rowDict)
    # close fiona object
    pointShp.close()
    print(f"----------- Vent file is saved as shape file in:'{shp_vent_file}'---------------")

def write_single_vent_shp(flow_id, long, lat, shp_vent_file, epsg_code):
    schema = {
        'geometry': 'Point',
        'properties': [("flow_id", 'str')],
    }

    with fiona.open(shp_vent_file, mode='w',
                    driver='ESRI Shapefile',
                    schema=schema,
                    crs=f'epsg:{epsg_code}') as shp:
        point_geom = Point(float(long), float(lat))
        record = {
            'geometry': mapping(point_geom),
            'properties': {'flow_id': flow_id}
        }
        shp.write(record)

def crop_asc_file(sim_asc, cropped_asc_file):
    with open(sim_asc) as file:
        header_lines = [next(file) for _ in range(6)]
        ncols = int(header_lines[0].split()[1])
        nrows = int(header_lines[1].split()[1])
        xllcorner = float(header_lines[2].split()[1])
        yllcorner = float(header_lines[3].split()[1])
        cellsize = float(header_lines[4].split()[1])
        nodata_value = float(header_lines[5].split()[1])
        # Read the values from the ASC file
        data_lines = [line.strip().split() for line in file]
        # Convert the data lines to a NumPy array
        data = np.array(data_lines, dtype=float)
    # Determine the index of the non-zeros lines and columns
    nonzero_rows, nonzero_cols = np.nonzero(data)
    # Calculate the limit of the crop
    min_row, max_row = np.min(nonzero_rows), np.max(nonzero_rows)
    min_col, max_col = np.min(nonzero_cols), np.max(nonzero_cols)
    # Defnie the values of the cropped data
    cropped_data = data[min_row:max_row + 1, min_col:max_col + 1]
    # Update the headers of the new cropped asc file
    cropped_nrows, cropped_ncols = cropped_data.shape
    cropped_xllcorner = xllcorner + min_col * cellsize
    cropped_yllcorner = yllcorner + (nrows - max_row - 1) * cellsize

    header_lines = [
        f"ncols {cropped_ncols}\n",
        f"nrows {cropped_nrows}\n",
        f"xllcorner {cropped_xllcorner}\n",
        f"yllcorner {cropped_yllcorner}\n",
        f"cellsize {cellsize}\n",
        f"nodata_value {nodata_value}\n"
    ]
    # write data in the new cropped asc file
    with open(cropped_asc_file, "w") as file:
        file.writelines(header_lines)
        for row in cropped_data:
            line = " ".join(str(value) if value is not None else str(nodata_value) for value in row)
            file.write(line + "\n")

def convert_to_tiff(cropped_asc_file, sim_tif_file):
    with rasterio.open(cropped_asc_file) as src:
        profile = src.profile.copy()
        profile["compress"] = "deflate"  # Use deflate compression
        profile["tiled"] = True  # Enable tiling for better performance and compression
        profile["blockxsize"] = 128  # Adjust the tile size as needed
        profile["blockysize"] = 128
        data = src.read(1)
        with rasterio.open(sim_tif_file, "w", **profile) as dst:
            dst.write(data, 1)

def crop_and_convert_to_tif(sim_asc, cropped_geotiff_file, epsg_code):
    """
    Crops an ASC file and saves it as a GeoTIFF file.

    :param sim_asc: Path to the input ASC file
    :param cropped_geotiff_file: Path to the output GeoTIFF file
    :param epsg_code: EPSG code for the coordinate reference system
    """
    # Read the ASCII file
    with open(sim_asc) as file:
        header_lines = [next(file) for _ in range(6)]
        ncols = int(header_lines[0].split()[1])
        nrows = int(header_lines[1].split()[1])
        xllcorner = float(header_lines[2].split()[1])
        yllcorner = float(header_lines[3].split()[1])
        cellsize = float(header_lines[4].split()[1])
        nodata_value = float(header_lines[5].split()[1])

        # Read the values from the ASC file
        data_lines = [line.strip().split() for line in file]
        # Convert the data lines to a NumPy array
        data = np.array(data_lines, dtype=float)

    # Determine the index of the non-zero rows and columns
    nonzero_rows, nonzero_cols = np.nonzero(data)
    # Calculate the limits of the crop
    min_row, max_row = np.min(nonzero_rows), np.max(nonzero_rows)
    min_col, max_col = np.min(nonzero_cols), np.max(nonzero_cols)

    # Define the values of the cropped data
    cropped_data = data[min_row:max_row + 1, min_col:max_col + 1]

    # Update the headers of the new cropped ASC file
    cropped_nrows, cropped_ncols = cropped_data.shape
    cropped_xllcorner = xllcorner + min_col * cellsize
    cropped_yllcorner = yllcorner + (nrows - max_row - 1) * cellsize

    # Define the transform and metadata for the GeoTIFF
    transform = from_origin(cropped_xllcorner, cropped_yllcorner + cropped_nrows * cellsize, cellsize, cellsize)
    metadata = {
        'driver': 'GTiff',
        'count': 1,
        'dtype': 'float32',
        'width': cropped_ncols,
        'height': cropped_nrows,
        'crs': CRS.from_epsg(epsg_code),
        'transform': transform,
        'nodata': nodata_value
    }


    # Write the data to a GeoTIFF file
    with rasterio.open(cropped_geotiff_file, 'w', **metadata) as dst:
        dst.write(cropped_data, 1)

    print(f" Cropped simulation saved in Geotiff at '{cropped_geotiff_file}'")

def get_runouts_grid_shp(run_outs_file, shp_runouts, epsg_code):

    pointDf = pd.read_csv(run_outs_file, header=0, sep=',')

    if pointDf.columns[0] == "":
        pointDf = pointDf.drop(pointDf.columns[0], axis=1)


    # 🔥 Detect effusion column automatically
    if "Effusion_rate" in pointDf.columns:
        eff_col = "Effusion_rate"
    elif "Effusion_r" in pointDf.columns:
        eff_col = "Effusion_r"
    elif "Eff_r" in pointDf.columns:
        eff_col = "Eff_r"
    elif "Effusion Rate" in pointDf.columns:
        eff_col = "Effusion Rate"
    else:
        raise ValueError(f"No effusion column found. Columns: {pointDf.columns}")

    # ✅ Mapping CSV → shapefile
    field_map = {
        eff_col: "Eff_r",
        "Median X": "Med_X",
        "Median Y": "Med_Y",
        "25th Percentile X": "P25_X",
        "25th Percentile Y": "P25_Y",
        "75th Percentile X": "P75_X",
        "75th Percentile Y": "P75_Y",
        "Average X": "Avg_X",
        "Average Y": "Avg_Y",
        "+2sigma X": "P2sig_X",
        "+2sigma Y": "P2sig_Y",
        "-2sigma X": "M2sig_X",
        "-2sigma Y": "M2sig_Y"
    }

    # ✅ Schema
    schema = {
        'geometry': 'Point',
        'properties': [(v, 'float') for v in field_map.values()]
    }

    # Force Eff_r as int
    schema['properties'][0] = ("Eff_r", 'int')

    with fiona.open(shp_runouts, mode='w',
                    driver='ESRI Shapefile',
                    schema=schema,
                    crs=f'epsg:{epsg_code}') as pointShp:

        for _, row in pointDf.iterrows():

            properties = {}

            for csv_col, shp_col in field_map.items():

                if csv_col not in row:
                    raise ValueError(f"Missing column: {csv_col}")

                value = row[csv_col]

                if csv_col == eff_col:
                    properties[shp_col] = int(value)
                else:
                    properties[shp_col] = float(value)

            rowDict = {
                'geometry': {
                    'type': 'Point',
                    'coordinates': (row["Median X"], row["Median Y"])
                },
                'properties': properties
            }

            pointShp.write(rowDict)

    print(f"----------- Runouts coordinates saved in: '{shp_runouts}' -----------")
def cut_lines_losd(losd_path, run_outs_path, output_path):
    """
    Cut LoSd lines into P25-P75 segments based on run_outs features,
    and save all segments in a single shapefile.

    Parameters:
        losd_path (str): Path to the LoSd line shapefile.
        run_outs_path (str): Path to run_outs shapefile with P25_X, P25_Y, P75_X, P75_Y, Eff_r.
        output_path (str): Path to save the output shapefile.

    Returns:
        None
    """
    # Load LoSd lines
    with fiona.open(losd_path, 'r') as losd_src:
        losd_lines = [shape(feat['geometry']) for feat in losd_src]

        crs = losd_src.crs  # keep the CRS
        schema = {
            'geometry': 'LineString',
            'properties': {'Eff_r': 'float'}
        }

    # Load run_outs features
    with fiona.open(run_outs_path, 'r') as runs_src:
        run_features = list(runs_src)

    all_segments = []

    for row in run_features:
        props = row['properties']
        eff_r = props['Eff_r']
        med_x = props.get('Med_X')
        med_y = props.get('Med_Y')
        p25 = Point(props['P25_X'], props['P25_Y'])
        p75 = Point(props['P75_X'], props['P75_Y'])


        for line in losd_lines:
            if not isinstance(line, LineString) or line.length == 0:
                continue

            # Project P25 and P75 onto the line
            d_start = line.project(p25)
            d_end = line.project(p75)

            if d_start > d_end:
                d_start, d_end = d_end, d_start

            # Extract segment
            coords = []
            for coord in line.coords:
                pt = Point(coord)
                d = line.project(pt)
                if d_start <= d <= d_end:
                    coords.append(coord)

            # Ensure segment starts/ends exactly at P25/P75
            coords = [line.interpolate(d_start).coords[0]] + coords + [line.interpolate(d_end).coords[0]]
            segment = LineString(coords)

            all_segments.append({
                'Eff_r': eff_r,
                'Med_X': med_x,
                'Med_Y': med_y,
                'geometry': segment
            })
    schema = {
        'geometry': 'LineString',
        'properties': {
            'Eff_r': 'float',
            'Med_X': 'float',
            'Med_Y': 'float'
        }
    }
    # Write output shapefile
    with fiona.open(output_path, 'w', driver='ESRI Shapefile', crs=crs, schema=schema) as dst:
        for seg in all_segments:
            dst.write({
                'geometry': mapping(seg['geometry']),
                'properties': {
                    'Eff_r': seg['Eff_r'],
                    'Med_X': seg['Med_X'],
                    'Med_Y': seg['Med_Y']
                }
            })

def cut_lines_losd_30pct(losd_path, run_outs_path, output_path):
    """
    Cut LoSd lines into segments centered on run_out points,
    with total length = 30% of distance_r.
    Save all segments in a single shapefile.

    Parameters:
        losd_path (str): Path to the LoSd line shapefile.
        run_outs_path (str): Path to run_outs shapefile with x_run_out, y_run_out, distance_r, Eff_r.
        output_path (str): Path to save the output shapefile.
    """
    # Load LoSd lines
    with fiona.open(losd_path, 'r') as losd_src:
        losd_lines = [shape(feat['geometry']) for feat in losd_src]
        crs = losd_src.crs

    # Load run_outs features
    with fiona.open(run_outs_path, 'r') as runs_src:
        run_features = list(runs_src)

    all_segments = []

    for row in run_features:
        props = row['properties']
        eff_r = props.get('Eff_r')
        x_run_out = props.get('X_run_out')
        y_run_out = props.get('Y_run_out')
        distance_r = props.get('Distance_r')

        if None in (eff_r, x_run_out, y_run_out, distance_r):
            print("⚠️ Skipping row, missing values:", props)
            continue

        midpoint = Point(x_run_out, y_run_out)
        seg_length = 0.3 * distance_r
        half_length = seg_length / 2.0

        for line in losd_lines:
            if not isinstance(line, LineString) or line.length == 0:
                continue

            # Project midpoint onto the line
            d_mid = line.project(midpoint)

            # Define start and end distances
            d_start = max(0, d_mid - half_length)
            d_end = min(line.length, d_mid + half_length)

            # Always interpolate start and end
            start_pt = line.interpolate(d_start).coords[0]
            end_pt = line.interpolate(d_end).coords[0]

            # Collect intermediate coords
            coords = [start_pt]
            for coord in line.coords:
                pt = Point(coord)
                d = line.project(pt)
                if d_start < d < d_end:
                    coords.append(coord)
            coords.append(end_pt)

            # Build the segment
            segment = LineString(coords)

            all_segments.append({
                'Eff_r': eff_r,
                'X_run_out': x_run_out,
                'Y_run_out': y_run_out,
                'Distance_r': distance_r,
                'geometry': segment
            })

    # Define schema
    schema = {
        'geometry': 'LineString',
        'properties': {
            'Eff_r': 'int',
            'X_run_out': 'float',
            'Y_run_out': 'float',
            'Distance_r': 'int'
        }
    }

    # Write output shapefile
    with fiona.open(output_path, 'w', driver='ESRI Shapefile', crs=crs, schema=schema) as dst:
        for seg in all_segments:
            dst.write({
                'geometry': mapping(seg['geometry']),
                'properties': {
                    'Eff_r': seg['Eff_r'],
                    'X_run_out': seg['X_run_out'],
                    'Y_run_out': seg['Y_run_out'],
                    'Distance_r': seg['Distance_r']
                }
            })

def get_vents_runouts_shp(output_csv, shp_vents_runouts, epsg_code):

    # import points from slope file
    with open(output_csv,'r') as file:
        try:
            dialect = csv.Sniffer().sniff(file.read(1024))
            file.seek(0)  # Reset file pointer after sniffing
        except csv.Error:
            print("Could not determine delimiter.")
            return

    df = pd.read_csv(output_csv, sep=dialect.delimiter)

    # Check required columns
    if 'X' not in df.columns or 'Y' not in df.columns:
        raise ValueError("CSV must contain 'X' and 'Y' columns")

    schema = {
        'geometry': 'Point',
        'properties': {}  # No attributes
    }

    with fiona.open(
        shp_vents_runouts,
        mode='w',
        driver='ESRI Shapefile',
        schema=schema,
        crs=f'epsg:{epsg_code}'
    ) as shp:

        for _, row in df.iterrows():
            shp.write({
                'geometry': {
                    'type': 'Point',
                    'coordinates': (float(row['X']), float(row['Y']))
                },
                'properties': {}
            })
    print(f"----------- Vent file from run_outs is saved as shape file in:'{shp_vents_runouts}'---------------")