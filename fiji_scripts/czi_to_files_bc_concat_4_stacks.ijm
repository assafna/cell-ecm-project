experiments_path = "G:/My Drive/BGU/Thesis/Cell-ECM & Cell-ECM-Cell Project/Data/Experiments/";
raw_path = experiments_path + "Raw/";
serieses_path = experiments_path + "Manipulations/Serieses/";
experiment = "SN26_BlebAdded";
first_file_path = raw_path + experiment + "/Rep_Loc_120to165minPostseed.czi";
second_file_path = raw_path + experiment + "/Rep_Loc_165minto4hPostseed.czi";
third_file_path = raw_path + experiment + "/Rep_Loc_4h30minPostseed_postBlebb.czi";
fourth_file_path = raw_path + experiment + "/Rep_Loc_4h30minto8hPostseed_postBlebb.czi";
experiment_serieses_path = serieses_path + experiment + "/";

for (s = 1; s <= 24; s++) {
	run("Bio-Formats", "open=[" + first_file_path + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + s);
	rename("first");
	run("Bio-Formats", "open=[" + second_file_path + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + s);
	rename("second");
	run("Bio-Formats", "open=[" + third_file_path + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + s);
	rename("third");
	
	run("Duplicate...", "duplicate frames=1");
	rename("third_new");
	
	run("Bio-Formats", "open=[" + fourth_file_path + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_" + s);
	rename("fourth");
	run("Concatenate...", "  title=new_image open image1=first image2=second image3=third_new image4=fourth image5=[-- None --]");
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