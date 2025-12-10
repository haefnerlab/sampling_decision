function [orientation_bins_center, orientation_energy] = get_ori_energy_model(gabor, num_ori_bin)
% GET_ORI_ENERGY_MODEL
%   Compute orientation energy (0..pi) of a 2D pattern using its 2D FFT.
%   - Removes DC
%   - Restricts to a radial frequency band (annulus)
%   - Normalizes by number of frequency samples per bin
%
%   Inputs:
%     gabor        : 2D image (e.g., Gabor)
%     num_ori_bin  : number of orientation bins between 0 and pi
%
%   Outputs:
%     orientation_bins_center : center of each orientation bin (radians, 0..pi)
%     orientation_energy      : orientation energy per bin (average power per freq sample)

% -------- 1. 2D FFT --------
F = fftshift(fft2(gabor));          % Compute 2D FFT and center it
magnitude = abs(F);                 % Amplitude spectrum

[Nr, Nc] = size(gabor);
if Nr ~= Nc
    error('Function currently assumes a square image (Nr == Nc).');
end
N = Nr;

% -------- 2. Frequency grid (aligned with fftshift) --------
% Frequencies go from -0.5 to 0.5 (Nyquist-normalized)
freq_axis = (-floor(N/2):ceil(N/2)-1) / N;   % length N
[fx, fy] = meshgrid(freq_axis, freq_axis);   % same size as FFT

% -------- 3. Orientation map (0..pi) --------
orientation = atan2(fy, fx);   % range [-pi, pi]
orientation = mod(orientation, pi);   % wrap to [0, pi)

% -------- 4. Remove DC (undefined orientation) --------
dc_idx = floor(N/2) + 1;
magnitude(dc_idx, dc_idx) = 0;

% -------- 5. Radial frequency mask (annulus) --------
r = sqrt(fx.^2 + fy.^2);   % normalized radius in frequency space

% Set these according to your stimulus:
% r_min > 0 to remove DC / very low freq
% r_max <= 0.5 (Nyquist); you can set it close to 0.5 or around your Gabor carrier
r_min = 0.05;
r_max = 0.45;

radial_mask = (r >= r_min) & (r <= r_max);

% Apply radial mask once to magnitude (so masked points are zero)
magnitude = magnitude .* radial_mask;

% -------- 6. Orientation bins (0..pi) --------
bin_edges = linspace(0, pi, num_ori_bin + 1) - pi/num_ori_bin/2;   % num_ori_bin bins, num_ori_bin+1 edges
%%%% make sure the center has 0, 45, 90, 135
orientation_bins_center = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

orientation_energy = zeros(1, num_ori_bin);

% -------- 7. Bin energy by orientation, normalize by count --------
mag2 = magnitude.^2;   % power/energy

for i = 1:num_ori_bin
    if i < num_ori_bin
        ori_mask = (orientation >= bin_edges(i)) & (orientation < bin_edges(i+1));
    else
        % Include the upper edge in the last bin
        ori_mask = (orientation >= bin_edges(i)) & (orientation <= bin_edges(i+1));
    end

    idx = ori_mask & radial_mask;

    % Sum energy in this bin
    e = sum(mag2(idx), 'all');

    % Normalize by number of points in this bin to reduce sampling bias
    n = nnz(idx);
    if n > 0
        orientation_energy(i) = e / n;
    else
        orientation_energy(i) = 0;
    end
    %orientation_energy(i)  = e;
end

end
