function out = BackwardsComp(out);
	%Manually handle variables whose names have changed
	out.Projection.nO = out.Projection.G.number_orientations;
	out.Projection.nL = out.Projection.G.number_locations;
	out.Projection.nX = out.Projection.G.dimension_X;
	out.Projection.nG = out.Projection.G.dimension_G;
	out.Projection.pT = out.Projection.G.prior_task;
	
	out.InputImage.dyn = out.Projection.I.stimulus_regime;
	out.InputImage.c = out.Projection.I.stimulus_contrast;
	out.Sampling.n_burn = out.Projection.S.number_burn_in;
	out.Sampling.n_use = out.Projection.S.number_samples_to_use;
	out.Sampling.nsamples_per_evidence = out.Projection.G.number_samples_per_evidence;
	out.Projection.delta = out.Projection.G.delta;
	out.Projection.kappa_O = out.Projection.G.kappa_O;
	out.Projection.kappa_G = out.Projection.G.kappa_G;


	%Automatically update non-name changing variables 
	names = fieldnames(out.Projection.G);
	data = struct2cell(out.Projection.G);
	l = length(names);
	for i=1:l;
		out.Projection.(names{i}) = data{i};
	end
	names = fieldnames(out.Projection.S);
	data = struct2cell(out.Projection.S);
	l = length(names);
	for i=1:l;
		out.Projection.(names{i}) = data{i};
	end
	names = fieldnames(out.Projection.I);
	data = struct2cell(out.Projection.I);
	l = length(names);
	for i=1:l;
		out.Projection.(names{i}) = data{i};
	end