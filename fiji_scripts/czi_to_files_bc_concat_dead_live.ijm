experiments_path = "G:/My Drive/BGU/Thesis/Cell-ECM & Cell-ECM-Cell Project/Data/Experiments/";
raw_path = experiments_path + "Raw/";
serieses_path = experiments_path + "Manipulations/Serieses/";
experiment = "LiveDead_201117";
first_file_path = raw_path + experiment + "/1/Live-dead Experiment_171120 (first time cycle).czi";
second_file_path = raw_path + experiment + "/2/Live-dead Experiment_171120.czi";
experiment_serieses_path = serieses_path + experiment + "/";

for (s = 1; s <= 13; s++) {
	run("Bio-Formats", "open=[" + first_file_path + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + s);
	rename("first");
	run("8-bit");
	run("Bio-Formats", "open=[" + second_file_path + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + s);
	rename("second");
	run("8-bit");
	run("Concatenate...", "  title=new_image open image1=first image2=second image3=[-- None --]");
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