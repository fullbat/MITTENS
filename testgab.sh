nifti_mitns = MITTENS(nifti_prefix="mgh")
loaded_voxel_graph = VoxelGraph("mgh_voxel_graph.mat")
paths, probs = voxel_graph.region_voxels_to_region_query("Frontal_Inf_Orb_R.nii.gz","Hippocampus_L.nii.gz", write_trk="test", write_prob="test")
