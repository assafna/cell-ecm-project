experiments_path = "G:/My Drive/BGU/Thesis/Cell-ECM & Cell-ECM-Cell Project/Data/Experiments/";
raw_path = experiments_path + "Raw/";
serieses_path = experiments_path + "Manipulations/Serieses/";
experiment = "SN20_Bleb_fromStart";
file_path = raw_path + experiment + "/SN20_rep3_Blebbistatin_from_start_time_laps_t1to8hh_tile.czi";
experiment_serieses_path = serieses_path + experiment + "/";

for (s = 1; s <= 20; s++) {
	run("Bio-Formats", "open=[" + file_path + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + s);
	rename("new_image");
	run("Duplicate...", "title=red duplicate channels=2");
	run("Grays");
	run("Bleach Correction", "correction=[Histogram Matching]");
	selectWindow("new_image");
	run("Duplicate...", "title=green duplicate channels=1");
	run("Green");
	run("Merge Channels...", "c1=DUP_red c2=green create");
	File.makeDirectory(experiment_serieses_path);
	saveAs("Tiff", experiment_serieses_path + "series_" + s + "_bc.tif");
	run("Close All");
}