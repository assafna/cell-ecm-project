experiment = "SN26_BlebAdded";

for (series = 1; series <= 24; series++) {
	open("G:/My Drive/BGU/Thesis/Cell-ECM & Cell-ECM-Cell Project/Data/Experiments/Manipulations/Serieses/" + experiment + "/series_" + series + "_bc.tif");
	experiment_objects_path = "G:/My Drive/BGU/Thesis/Cell-ECM & Cell-ECM-Cell Project/Data/Experiments/Manipulations/Objects/" + experiment;
	File.makeDirectory(experiment_objects_path);
	series_objects_dir = experiment_objects_path + "/Series " + series;
	File.makeDirectory(series_objects_dir);

	info = getImageInfo();
	lines = split(info, "\n");
	for (line_index = 0; line_index < lines.length; line_index++) {
		line = lines[line_index];
		if (startsWith(line, "  Frame: ")) {
			frames = split(line, "/");
			tps = parseInt(frames[1]);
		}
	}
	
	for (tp = 1; tp <= tps; tp++) {
		run("Duplicate...", "duplicate channels=2 frames=" + tp);

		// check if the image is black
		run("Set Measurements...", "mean redirect=None decimal=3");
		run("Measure");
		mean = getResult("Mean");
		run("Close");
		if (mean == 0) {
			File.saveString("BLACK", series_objects_dir + "/tp_" + tp + ".txt");
		} else {
			run("Median...", "radius=2 stack");
			run("Gaussian Blur...", "sigma=2 scaled stack");
			run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=10 show_numbers white_numbers store_results_within_a_table_named_after_the_image_(macro_friendly) redirect_to=none");
			run("3D Objects Counter", "threshold=15 slice=18 min.=400 max.=9437184 objects statistics");
			saveAs("Results", series_objects_dir + "/tp_" + tp + ".txt");
			if (tp == 1) {
				run("glasbey inverted");
				saveAs("Tiff", series_objects_dir + "/objects_map.tif");	
			}	
			selectWindow("tp_" + tp + ".txt");
			run("Close");
			close();
			close();
		}
		selectWindow("series_" + series + "_bc.tif");
	}
	
	run("Close All");	
}