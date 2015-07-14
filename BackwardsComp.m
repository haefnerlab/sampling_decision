function out = BackwardsComp(out);
	out.Projection.nO = out.Projection.G.number_orientations;
	out.Projection.nL = out.Projection.G.number_locations;
	out.Projection.nX = out.Projection.G.dimension_X;
	out.Projection.nG = out.Projection.G.dimension_G;
	out.Projection.pT = out.Projection.G.prior_task;
	
	out.InputImage.dyn = out.Projection.I.stimulus_regime;
	out.InputImage.c = out.Projection.I.stimulus_contrast;
	out.Sampling.n_burn = out.Projection.S.number_burn_in;
	out.Sampling.n_use = out.Projection.S.number_samples_to_use;
	out.Sampling.nsamples_per_evidence = out.Projection.S.number_samples_per_evidence;
