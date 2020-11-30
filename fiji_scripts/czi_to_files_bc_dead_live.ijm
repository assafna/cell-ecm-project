experiments_path = "G:/My Drive/BGU/Thesis/Cell-ECM & Cell-ECM-Cell Project/Data/Experiments/";
raw_path = experiments_path + "Raw/";
serieses_path = experiments_path + "Manipulations/Serieses/";
experiment = "LiveDead_201115";
file_path = raw_path + experiment + "/Live-dead Experiment_151120 (finished at 1800).czi";
experiment_serieses_path = serieses_path + experiment + "/";

for (s = 1; s <= 17; s++) {
	run("Bio-Formats", "open=[" + file_path + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + s);
	rename("new_image");
	run("8-bit");
	run("Duplicate...", "title=red duplicate channels=1");
	run("Grays");
	run("Bleach Correction", "correction=[Histogram Matching]");
	selectWindow("new_image");
	run("Duplicate...", "title=live duplicate channels=2");
	selectWindow("new_image");
	run("Duplicate...", "title=dead duplicate channels=3");
	imageCalculator("Add create stack", "live","dead");
	selectWindow("Result of live");
	rename("green");
	run("Green");
	run("Merge Channels...", "c1=DUP_red c2=green create");
	File.makeDirectory(experiment_serieses_path);
	saveAs("Tiff", experiment_serieses_path + "series_" + s + "_bc.tif");
	run("Close All");
}