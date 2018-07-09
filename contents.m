% BSIE: Blind System Identification and Equalization Toolbox
%
% Initialize Blind Adaptive System Identification
%   init_mclms      - Initialize Multichannel Least Mean Squares
%   init_mcn        - Initialize Multichannel Newton
%   init_mcflms     - Initialize Multichannel Frequency Domain LMS
%   init_nmcflms    - Initialize Normalized MCFLMS
%   init_nmcflms_sc - Initialize Normalized MCFLMS with Spectral Constraint
%   init_nmcmdflms  - Initialize Normalized Multichannel MDF LMS
%   init_rnmcflms   - Initialize Robust NMCFLMS
%
% Blind System Identification
%   mcss            - Multichannel Subspace Algorithm
%
% Blind Adaptive System Identification
%   mclms           - Multichannel Least Mean Squares Algorithm
%   mcn             - Multichannel Newton Algorithm
%   mcflms          - Multichannel Frequency Domain LMS Algorithm
%   nmcflms         - Normalized MCFLMS Algorithm
%   nmcflms_sc      - Normalized MCFLMS with Spectral Constraint Algorithm
%   nmcmdflms       - Normalized Multichannel MDF LMS Algorithm
%   rnmcflms        - Robust NMCFLMS Algorithm
%   
% System Equalization
%   lsinvfilt           - Compute Least Squares Inverse Filter
%   wls                 - Compute Weighted LS Filter
%   wls_iterative       - Compute Iterative Weighted LS Filter
%   channel_shortening  - Compute Channel Shortening Filter
%   rmcls               - Compute Relaxed Multichannel LS Filter
%   rcic                - Compute RMCLS Filter with Constrained Initial Coefficients
%   rec                 - Compute RMCLS Filter with Envelope Constraint
%
% Performance Evaluation and Analyses
%   edc                 - Compute Energy Decay Curve
%   npm                 - Normalized Projection Misalignment
%   npm_ac              - NPM with Alignment Correction
%   magnitude_deviation - Compute Magnitude Deviation
%   phase_deviation     - Compute Phase Distortion
%   gmc_st              - Generalized Multichannel Clustering
%
% Auxiliary Functions
%   generate_data       - Generate Sensor Data
%   generate_sie        - Generate System Identification Errors
%   rdp                 - Remove Direct-Path Propagation
%   ula_pos             - Sensor Positions of a Uniform Linear Array
%
%   Copyright (C) Imperial College London 2009-2011
%
%   BSIE is a MATLAB toolbox for blind SIMO identification and equalization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   ftp://prep.ai.mit.edu/pub/gnu/COPYING-2.0 or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%