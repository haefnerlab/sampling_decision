function out = BackwardsComp(out);
	out.Projection.nO = out.Projection.number_orientations;
	out.Projection.nL = out.Projection.number_locations;
	out.Projection.nX = out.Projection.dimension_X;
	out.Projection.nG = out.Projection.dimension_G;
	out.Projection.pT = out.Projection.prior_task;
	
	out.InputImage.dyn = out.Projection.stimulus_regime;
	out.InputImage.c = out.Projection.stimulus_contrast;
	out.Sampling.n_burn = out.Projection.number_burn_in;
	out.Sampling.n_use = out.Projection.number_samples_to_use;
	out.Sampling.nsamples_per_evidence = out.Projection.number_samples_per_evidence;
