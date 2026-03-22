import fiona
import matplotlib
import matplotlib.pyplot as plt
from shapely.geometry import shape, Point
import matplotlib.colors as colors
import rasterio
from PIL import Image
from adjustText import adjust_text
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
import os
import datetime
import numpy as np
from rasterio.transform import Affine
from rasterio.crs import CRS
from matplotlib.lines import Line2D
import matplotlib.cm as cm

matplotlib.rcParams['figure.dpi'] = 150        # for display in interactive window
matplotlib.rcParams['savefig.dpi'] = 300       # for saving with the GUI or savefig()

class Mapping:
    def __init__(self, path_to_folder, dem, flow_id, map_layers, sim_layers, mode, language, grid_mode):
        """
        Create a map for either 'downflow' or 'downflowgo' based on the mode provided and in English or in French

        Parameters:
        - path_to_folder: Folder to save the generated map.
        - dem: DEM file path.
        - flow_id: ID of the flow.
        - map_layers: Dictionary containing map layer file paths.
        - sim_layers: Dictionary containing simulation layer file paths.
        - mode: Either 'downflow' or 'downflowgo' to choose the type of map. Default is 'downflow'.
        """
        plt.close('all')

        self.path_to_folder = path_to_folder
        self.dem = dem
        self.flow_id = flow_id
        self.mode = mode
        self.language = language
        self.grid_mode = grid_mode

        # Create the map figure
        self.fig, self.ax = plt.subplots(figsize=(9, 5))

        current_datetime = datetime.datetime.now()
        dt = datetime.datetime.now(datetime.timezone.utc)
        Date = dt.strftime('%d-%m-%y %H:%M UTC')
        Date_str_fr = dt.strftime('%d-%m-%y à %H:%M UTC')
        Date_str_en = dt.strftime('%d-%m-%y at %H:%M UTC')


        self.losd_path = sim_layers['shp_losd_file']
        self.vent_path = sim_layers['shp_vent_file']
        self.sim_tif_file = sim_layers['cropped_geotiff_file']
        self.shp_iqr = None
        self.shp_30pct = None

        if mode == 'downflowgo':
            self.run_outs_path = sim_layers['shp_runouts']
            if grid_mode == 'yes':
                self.shp_iqr = sim_layers['shp_iqr']
            if grid_mode == 'no':
                self.shp_30pct = sim_layers.get('shp_30pct')

        self.tiff_file = map_layers['img_tif_map_background']
        self.source_img_tif_map_background = map_layers["source_img_tif_map_background"]
        self.lavaflow_outline_path = map_layers['lava_flow_outline_path']
        self.monitoring_network_path = map_layers['monitoring_network_path']
        self.logo = map_layers['logo_path']
        self.unverified_data = map_layers['unverified_data']

        # Label part
        self.probable_events_label = {
            'FR':'Events',
            'EN':'Vents'
        }

        self.main_vent_label = {
            'FR':'Bouche éruptive ('+ flow_id+')',
            'EN':'Main vent (' + flow_id + ')'
        }

        self.steepest_descend_line_label = {
            'FR':'Trajectoire principale',
            'EN':'Line of steepest descent'
        }

        self.median_runouts_label = {
            'FR':'Distances atteintes \npour un débit donné (m$^{3}$/s),\n intervalle interquartile',
            'EN':'Median runouts with interquartiles range \nfor a given effusion rate (m$^{3}$/s)'
        }

        self.runouts_label = {
            'FR':'Distances atteintes (\u00B1 30%) \npour un débit donné (m$^{3}$/s)',
            'EN':'Runouts (\u00B1 30%) for a given \neffusion rate (m$^{3}$/s)'
        }

        self.monitoring_network_label = {
            'FR':'Stations de surveillance',
            'EN':'Monitoring network'
        }

        self.lava_flow_outline_label = {
            'FR':'Contour de la coulée',
            'EN':'Lava flow outline'
        }

        self.title = {
            'FR':f"Carte de simulation {'DOWNFLOWGO' if mode=='downflowgo' else 'DOWNFLOW'} pour l\'éruption en cours",# + str(Date),
            'EN':f"Short-term hazard map using {'DOWNFLOWGO' if mode=='downflowgo' else 'DOWNFLOW'} modelling \n for the ongoing eruption"# + str(Date)
        }
        
        self.lava_flow_proba_label = {
            'FR':"Probabilité de passage de la coulée",
            'EN':"Lava flow path probability"
        }

        self.low_label = {
            'FR':"Basse",
            'EN':"Low"
        }

        self.high_label = {
            'FR':"Haute",
            'EN':"High"
        }
        self.map_date = {
            'FR':"Carte réalisée le "+ Date_str_fr,
            'EN':"Map made on "+ Date_str_en
        }

    def plot_background(self):
        # Load the TIFF image, convert it to RGB and read
        if self.tiff_file and self.tiff_file != "0":
            tiff_image = Image.open(self.tiff_file)
            tiff_image_rgb = tiff_image.convert('RGB')
            with rasterio.open(self.tiff_file) as src:
                tiff_image = src.read(1)
                extent = src.bounds
            # Display the TIFF image as the background
            tif = self.ax.imshow(tiff_image_rgb, extent=(extent.left, extent.right, extent.bottom, extent.top))

    def plot_simulation(self):
        # Load the compressed GeoTIFF file and extract the necessary information and data
        with rasterio.open(self.sim_tif_file) as src:
            data = src.read(1)
            transform = src.transform

        max_prob = int(data.max())
        # Define value intervals and associated colors
        if self.grid_mode == 'yes':
            data_min = 1e-9
            data_max = 0.0005
            cmap_base = cm.get_cmap('YlOrRd')
            colors_with_zero = [(0, 0, 0, 0)] + [cmap_base(i) for i in range(cmap_base.N)]
            cmap = colors.ListedColormap(colors_with_zero)
            self.intervals = np.linspace(data_min, data_max, cmap_base.N)
            self.intervals = np.insert(self.intervals, 0, 0)
            self.intervals = np.append(self.intervals, np.inf)
            norm = colors.Normalize(vmin=data_min, vmax=data_max, clip=True)

        else:
            # Cas max_prob > 1 : puissances de 10
            max_power = int(np.ceil(np.log10(max_prob)))
            self.intervals = [0] + [10 ** i for i in range(0, max_power + 1)]
            n_bins = len(self.intervals) - 1

            cmap_base = cm.get_cmap('YlOrRd')
            colors_with_zero = [(0, 0, 0, 0)] + [cmap_base(i) for i in range(cmap_base.N)]
            cmap = colors.ListedColormap(colors_with_zero)

            #gradient_colors = plt.cm.Reds(np.linspace(0.4, 0.9, n_bins))
            #gradient_colors[0] = (1, 0.84, 0, 1)  # gold
            #if n_bins > 1:
            #    gradient_colors[1] = (1, 0.65, 0, 1)  # orange

            #colors_list = [(0, 0, 0, 0)] + [tuple(c) for c in gradient_colors]
            #cmap = colors.ListedColormap(colors_list, name="my_raster_colormap")
            norm = colors.BoundaryNorm(self.intervals, cmap.N)

        # Plot the simulation data on the map
        self.ax.imshow(data, extent=(transform[2], transform[2] + transform[0] * data.shape[1],
                                    transform[5] + transform[4] * data.shape[0], transform[5]),  # Flip the y-axis
                        cmap=cmap, norm=norm)

        # Set the limit of the figure based on used DEM
        with open(self.dem) as file:
            header_lines = [next(file) for _ in range(6)]
            ncols_dem = int(header_lines[0].split()[1])
            nrows_dem = int(header_lines[1].split()[1])
            xllcorner_dem = float(header_lines[2].split()[1])
            yllcorner_dem = float(header_lines[3].split()[1])
            cellsize_dem = float(header_lines[4].split()[1])
            nodata_value = float(header_lines[5].split()[1])

        self.ax.set_xlim(xllcorner_dem, xllcorner_dem + ncols_dem * cellsize_dem)
        self.ax.set_ylim(yllcorner_dem, yllcorner_dem + nrows_dem * cellsize_dem)

    def plot_vector_layers(self):
        # Plot the vent vector layer
        with fiona.open(self.vent_path) as vent:
            for feature in vent:
                geometry = shape(feature['geometry'])
                if geometry.geom_type == 'Point':
                    x, y = geometry.x, geometry.y
                    if self.grid_mode=='yes':
                        self.ax.plot(x, y, 'yo', markersize=1, zorder=1)
                    else:
                        self.ax.plot(x, y, 'r^', markersize=10,zorder=1)

        # Plot the LOSD vector layer
        with fiona.open(self.losd_path) as losd:
            for feature in losd:
                geometry = shape(feature['geometry'])
                if geometry.geom_type == 'LineString':
                    x, y = zip(*geometry.coords)
                    self.ax.plot(x, y, color='#8B0000', linestyle='-', linewidth=2)


        ##Plot additional interquartil rectangular for grid:
        #if grid_mode == 'yes':
        #    tolerance = 50  # distance minimale in meters to avoid overlap
        #    eff_r_values = []
        #    with fiona.open(shp_iqr) as iqr:
        #        for feature in iqr:
        #            eff_r_values.append(feature['properties']['Effusion_r'])
        #
        #    # Normalisation for colormap
        #    norm_eff_r = colors.Normalize(vmin=min(eff_r_values), vmax=max(eff_r_values))
        #    #≈cmap_eff_r = cm.get_cmap('Blues') # change colormap if needed
        #    cmap_eff_r = lambda x: 'cyan'  # returns 'blue' for any input
        #
        #    seen_lines = []
        #    with fiona.open(shp_iqr) as iqr:
        #        for feature in iqr:
        #            geometry = shape(feature['geometry'])
        #            props = feature['properties']
        #            eff_r = props['Effusion_r']
        #            color = cmap_eff_r(norm_eff_r(eff_r))
        #
        #            if geometry.geom_type == 'LineString':
        #                # Avoid tracing two times the same
        #                is_duplicate = any(geometry.distance(s) < tolerance for s in seen_lines)
        #                if not is_duplicate:
        #                    x, y = zip(*geometry.coords)
        #                    ax.plot(x, y, color=color, linewidth=5, alpha=0.8) # plot the lines for a given effusion rate
        #                    seen_lines.append(geometry)
        #
        #                    # 2) plot Med_X / Med_Y and draw the arrow with the label
        #                    if props.get('Med_X') is not None and props.get('Med_Y') is not None:
        #                        mx, my = props['Med_X'], props['Med_Y']
        #                        ax.plot(mx, my,  'b', marker=7)
        #                        #ax.plot(mx, my, 'b', marker=7)
        #
        #                        # Add annotation with line and white box around label
        #                        ax.annotate(
        #                            str(int(eff_r)),
        #                            xy=(mx, my),  # tip of the arrow
        #                            xytext=(mx-200, my+800),  # base of the arrow where the label is
        #                            arrowprops=dict(arrowstyle='-',linewidth=1, color='blue'
        #                            ),
        #                            fontsize=8,
        #                            color='blue',
        #                            ha='center',
        #                            va="center",
        #                            bbox=dict(
        #                                facecolor='white',
        #                                edgecolor='blue',
        #                                alpha=0.9,
        #                                boxstyle='square,pad=0.2'
        #                            )
        #                        )
        #
        #    # Plot additional runout layers for 'downflowgo'
        #if grid_mode !='yes' and mode == 'downflowgo' and run_outs_path:
        #    with fiona.open(run_outs_path) as run_outs:
        #        points = []
        #        labels = []
        #
        #        for feature in run_outs:
        #            geometry = shape(feature['geometry'])
        #            properties = feature['properties']
        #            if geometry.geom_type == 'Point':
        #                x, y = geometry.x, geometry.y
        #                points.append((x, y))
        #                labels.append(properties['Effusion_r'])
        #                ax.plot(x, y, 'b', marker=7)
        #
        #        unique_points = {}
        #        for (x, y), label in zip(points, labels):
        #            if (x, y) not in unique_points or label < unique_points[(x, y)]:
        #                unique_points[(x, y)] = label
        #
        #        for (x, y), label in unique_points.items():
        #            ax.annotate(label, (x, y), xytext=(0, 8), textcoords='offset points', color='blue', weight='bold',
        #                        fontsize=10, ha='center')

        tolerance = 50  # distance minimale in meters to avoid overlap
        eff_r_values = []

        # Select shapefile depending on grid_mode
        shp_to_open = None
        if self.mode == 'downflowgo':
            if self.grid_mode == 'no':
                shp_to_open = self.shp_30pct

        if self.grid_mode == 'yes':
            shp_to_open = self.shp_iqr

        if shp_to_open is not None:
            with fiona.open(shp_to_open) as src:
                for feature in src:
                    eff_r_values.append(feature['properties']['Eff_r'])

            # Normalisation for colormap
            norm_eff_r = colors.Normalize(vmin=min(eff_r_values), vmax=max(eff_r_values))
            # cmap_eff_r = cm.get_cmap('Blues') # change colormap if needed
            cmap_eff_r = lambda x: 'cyan'

            seen_lines = []

            with fiona.open(shp_to_open) as src:
                for feature in src:
                    geometry = shape(feature['geometry'])
                    props = feature['properties']
                    eff_r = props['Eff_r']
                    color = cmap_eff_r(norm_eff_r(eff_r))

                    if geometry.geom_type == 'LineString':
                        # Avoid tracing two times the same
                        is_duplicate = any(geometry.distance(s) < tolerance for s in seen_lines)
                        if not is_duplicate:
                            x, y = zip(*geometry.coords)
                            self.ax.plot(x, y, color=color, linewidth=5, alpha=0.8)
                            seen_lines.append(geometry)

                            if self.grid_mode == 'yes':
                                mx, my = props.get('Med_X'), props.get('Med_Y')
                            else:
                                mx, my = props.get('X_run_out'), props.get('Y_run_out')

                            if mx is not None and my is not None:
                                self.ax.plot(mx, my, 'b', marker=7)

                                # Add annotation with line and white box around label
                                self.ax.annotate(
                                    str(int(eff_r)),
                                    xy=(mx, my),  # tip of the arrow
                                    xytext=(mx - 200, my + 800),  # base of the arrow
                                    arrowprops=dict(arrowstyle='-', linewidth=1, color='blue'),
                                    fontsize=8,
                                    color='blue',
                                    ha='center',
                                    va="center",
                                    bbox=dict(
                                        facecolor='white',
                                        edgecolor='blue',
                                        alpha=0.9,
                                        boxstyle='square,pad=0.2'
                                    )
                                )
        else:
            # No flowgo shapefile in this mode
            pass
            # ------------ Plot lava flow outline and monitoring network ----------

        if self.lavaflow_outline_path and self.lavaflow_outline_path != "0":
            with fiona.open(self.lavaflow_outline_path) as lavaflow_outline:
                for feature in lavaflow_outline:
                    if feature.get('geometry') is not None:
                        geometry = shape(feature['geometry'])
                        if geometry.geom_type == 'Polygon':
                            coordinates = geometry.exterior.coords.xy
                            x = coordinates[0]
                            y = coordinates[1]
                            polygon_coords = list(zip(x, y))
                            polygon = Polygon(polygon_coords, edgecolor='black', facecolor='none')
                            self.ax.add_patch(polygon)
                    else:
                        print("Lavaflow_outline does not have a valid geometry")

        # Plot monitoring network vector layer
        if self.monitoring_network_path and self.monitoring_network_path != "0":
            with fiona.open(self.monitoring_network_path) as monitoring_network:
                label_monitoring_network = []
                point_monitoring_network = []
                for feature in monitoring_network:
                    geometry = shape(feature['geometry'])
                    properties = feature['properties']
                    if geometry.geom_type == 'Point':
                        x, y = geometry.x, geometry.y
                        point_monitoring_network.append((x, y))
                        label_monitoring_network.append(properties['Name'])
                        self.ax.plot(x, y, color="#333333", linestyle='', marker='s', markersize=2)
                    else:
                        print("Monitoring_network does not have a valid geometry")
                # Create a dictionary to keep only the smallest label (lexicographically) for each unique point
                unique_monitoring_network = {}
                for (x, y), label in zip(point_monitoring_network, label_monitoring_network):
                    if label is None:
                        continue  # Skip stations with no name
                    if (x, y) not in unique_monitoring_network or label < unique_monitoring_network[(x, y)]:
                        unique_monitoring_network[(x, y)] = label

                # Displaying the monitoring_network names just above the points
                for (x, y), label in unique_monitoring_network.items():
                    self.ax.annotate(
                        label, (x, y),
                        xytext=(0, 3),  # Offset of 5 units above the point
                        textcoords='offset points',
                        color="#333333", fontsize=7,
                        ha='center'  # Center the text horizontally
                    )

    def add_legend(self):
        # Adjust the figure size
        self.fig.subplots_adjust(right=0.7)
        self.fig.subplots_adjust(left=0.1)
        self.fig.subplots_adjust(top=0.9)
        self.fig.subplots_adjust(bottom=0.1)

        # Add the legend image to the map
        if self.grid_mode == 'yes':
            self.ax.plot([], [], 'yo', markersize=3, label=self.probable_events_label[self.language])

        else:
            self.ax.plot([], [], 'r^', markersize=7, label=self.main_vent_label[self.language])
        self.ax.plot([], [], color='#8B0000', linestyle='-', linewidth=2, label=self.steepest_descend_line_label[self.language])

        if self.mode == 'downflowgo' and self.grid_mode == 'yes': #marker=r'$\downarrow$'
            #ax.plot([], [], 'bv', linestyle='None', markersize=5, label='Distances atteintes \npour un débit donné (m$^{3}$/s)')
            #ax.plot([], [], '-',color='cyan', linewidth=5, label='Variabilité médiane de la distance atteinte')
            self.ax.plot([], [],  marker='v', markersize=7,
                            markerfacecolor='blue', markeredgecolor='none',
                    color='cyan', linewidth=5,
                    label=self.median_runouts_label[self.language])
        if self.mode == 'downflowgo' and self.grid_mode == 'no':
            self.ax.plot([], [],  marker='v', markersize=7,
                            markerfacecolor='blue', markeredgecolor='none',
                    color='cyan', linewidth=5,
                    label=self.runouts_label[self.language])
            #ax.plot([], [], '-',color='cyan', linewidth=5, marker='bv', markersize=5, label=' 30 % incertitude sur la distance')

        # Plot monitoring network vector layer
        if self.monitoring_network_path and self.monitoring_network_path != "0":
            #self.ax.plot([], [], color="#333333", linestyle='', marker='s', markersize=2, label=self.monitoring_network_label[self.language])
            self.ax.plot([], [], 'ks', markersize=2,label=self.monitoring_network_label[self.language])
        if self.lavaflow_outline_path and self.lavaflow_outline_path != "0":
            self.ax.plot([], [], 'k-', markersize=10, label=self.lava_flow_outline_label[self.language])

    def final_adjusment(self):
        #show label from yaxis en completo
        #formatter = ScalarFormatter(useMathText=False)
        #formatter.set_scientific(False)
        #ax.yaxis.set_major_formatter(formatter)

        # Add title and date
        self.ax.set_title(self.title[self.language])

        # Load Logo image
        if self.logo and self.logo != "0":
            img = plt.imread(self.logo)
            logo_box = OffsetImage(img, zoom=0.05)  # size of the logo
            logo_anchor = (1.3, 0.02)  # Define the coordinate of the image anchor
            # Create the annotation of the image (remove frame
            ab = AnnotationBbox(logo_box, logo_anchor, xycoords='axes fraction', frameon=False)
            self.ax.add_artist(ab)
        legend = self.ax.legend(bbox_to_anchor=(1, 0.7), loc="upper left", fontsize='8')
        legend.get_frame().set_linewidth(0)
        self.fig.canvas.draw()  # needed so legend position is finalized
        legend_box = legend.get_window_extent(self.fig.canvas.get_renderer())
        legend_box = legend_box.transformed(self.fig.transFigure.inverted())  # convert to figure coords
        
        # Create a colorbar for the simulation data in raster
        #cax = fig.add_axes([0.71, 0.70, 0.2, 0.02])  # Coordinates and size of the axis
        # Create the colorbar with the colors of the intervals

        # --- Colormap for colorbar (continuous fade gold→darkred) ---
        cmap_colorbar = cm.get_cmap('YlOrRd')
        norm2 = colors.Normalize(vmin=min(self.intervals), vmax=max(self.intervals))
        # Use LogNorm so colorbar matches log-like spacing
        #norm2 = colors.LogNorm(vmin=1, vmax=max(intervals))
        sm = ScalarMappable(cmap=cmap_colorbar, norm=norm2)
        sm.set_array([])  # Set empty array for the colorbar
        # Rotate label of Y and orientate vertically
        self.ax.set_yticks(self.ax.get_yticks())
        self.ax.set_yticklabels(self.ax.get_yticklabels(), rotation='vertical')
        cax = self.fig.add_axes([
            legend_box.x0,  # same left as legend
            legend_box.y1 + 0.05,  # just above legend (0.02 = small gap)
            0.2,  0.02  # width and height of colorbar
        ])
        cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')

        cbar.ax.set_title(self.lava_flow_proba_label[self.language], fontsize=8, loc='left', pad=14)
        cbar.ax.text(0, 1.5, self.low_label[self.language], fontsize=8, ha="left", va="center", transform=cbar.ax.transAxes)
        cbar.ax.text(1, 1.5, self.high_label[self.language], fontsize=8, ha="right", va="center", transform=cbar.ax.transAxes)
        cbar.ax.set_xticks([])  # clear ticks
        cbar.ax.set_xticklabels([])  # clear labels
        
        # Adding the "UNVERIFIED DATA" text in the middle of the map
        if self.unverified_data != '0':
            self.ax.text(
                0.5, 0.5,  self.unverified_data,#.replace("-", "\n"),
                transform=self.ax.transAxes,  # This makes the text relative to the axes
                fontsize=30, color='red', alpha=0.5,  # Red color, semi-transparent
                ha='center', va='center',  # Centered horizontally and vertically
                rotation=45,  # Rotate the text 45 degrees
                fontweight='bold')  # Bold text
            
        # Adding text for credit of the background map
        if self.source_img_tif_map_background != '0':
            self.ax.text(
                0.01, 0.01, self.source_img_tif_map_background,
                transform=self.ax.transAxes,
                fontsize=8, color='k',
                ha='left', va='bottom')

        # Adding text for map_date
        self.ax.text(1.46, 0.1,# Define the coordinate of the image anchor
                self.map_date[self.language],
                transform=self.ax.transAxes,
                fontsize=8, color='k',
                ha='right', va='bottom')

        # Add North arrow
        self.ax.annotate(
            "N",
            xy=(0.05, 0.95),  # arrow tip (in axes fraction coords)
            xytext=(0.05, 0.85),  # text position (below tip)
            xycoords="axes fraction",  # relative to axes (0–1)
            textcoords="axes fraction",
            ha="center", va="center",
            fontsize=12, fontweight="bold",
            arrowprops=dict(facecolor="black", headwidth=8)
        )
    
    def savepath(self):
        final_map = plt.savefig(os.path.join(self.path_to_folder, f"map_{self.flow_id}.png"), dpi=300, bbox_inches='tight')
        print('Map created',os.path.join(self.path_to_folder, f"map_{self.flow_id}.png"))
        return final_map
    
    def show_map(self, display):
        # Show if requested
        if display == 'yes':
            plt.show()

        # Close it if shown (good for scripts), or let user manage it if not
        if display == 'yes':
            plt.close(self.fig)


    def create_map(self, display):

        self.plot_background()

        # ------------ plot simulation .tiff ----------

        self.plot_simulation()

        # ------------ Plot vector layers simulation path + vent + run out ----------

        self.plot_vector_layers()

        # ------------ Add legend and colorbar ------------

        self.add_legend()

        # ------------ Final adjustments ------------

        self.final_adjusment()
        
        final_map = self.savepath()

        self.show_map(display)

        # Return the figure for use elsewhere (e.g., in notebook)
        return final_map