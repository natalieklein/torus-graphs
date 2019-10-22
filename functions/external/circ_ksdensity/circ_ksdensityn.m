function [vfEstimate, vfBinVol] = circ_ksdensityn(mfObservations, mfPDFSamples, mfDomains, vfSigmas, vfWeights)

% circ_ksdensityn FUNCTION - Compute a kernel density estimate over periodic and aperiodic domains
%
% Usage: [vfEstimate, vfBinVol] = circ_ksdensityn(mfObservations, mfPDFSamples, 
%                                                 <mfDomains, vfSigmas, vfWeights>)
%
% This function calculates a kernel density estimate of an (optionally
% weighted) data sample, over periodic and aperiodic domains. The sample is
% assumed to be independent across dimensions; i.e. density estimation is
% performed independently for each dimension of the data.
%
% 'mfObservations' is a set of observations made over a (possibly periodic)
% domain. Each row corresponds to a single observation, each column
% corresponds to a particular dimension. By default all dimensions are
% periodic in [0..2*pi]; this can be modified by providing the optional
% arument 'mfDomains'. Each row in 'mfDomains' is [fMin fMax], one row for
% each dimension in 'mfObservations'. If a particular dimension should not
% be periodic, the corresponding row should be [nan nan]. Bounded support
% over a dimension is NOT implemented; each dimension is either linear and
% infinite or periodic.
%
% 'mfPDFSamples' defines the sample points over which to perform
% the kernel density estimate, over the same domains as 'mfObservations'.
%
% Weighted estimations can be performed by providing the optional argument
% 'vfWeights', where each element in 'vfWeights' corresponds to the
% matching observation in 'mfObservations'.
%
% The kernel density estimate will be performed using a multivariate Gaussian
% kernel, independent along each dimension, and wrapped along the periodic
% dimensions as appropriate.  Kenel widths over periodic dimensions are estimated as
%    (4/3)^0.2 * circ_std(mfObservations(:, nDim), vfWeights) * (length(mfObservations)^-0.2)
%
% Kernel widths over non-periodic dimensions are estimated as
%    (4 * std(mfObservations(:, nDim), vfWeights)^5 / 3 / length(mfObservations))^(1/5)
%
% The optional argument 'fSigma' can be provided to set the width of the
% kernel.
%
% 'vfEstimate' will be a vector with a (weighted) histogram estimate of the
% underlying distribution, with an entry for each point in 'mfPDFSamples'.
% If no weighting is supplied, the estimate will be scaled to estimate a
% PDF over the supplied multi-dimensional domain, taking into account the
% estimated volume of each bin. If a weight vector is supplied, the
% estimate will be scaled such that the sum over the domain attempts to
% match the sum of weights, taking into account the multi-dimensional bin
% volumes.
%
% 'vfBinVol' is a vector containing volume estimates for each row in
% 'mfPDFSamples', under the assumption that each dimension is linearly
% scaled and mutually orthogonal.

% Author: Dylan Muir <dylan.muir@unibas.ch>
% Created: 23rd October, 2013

% -- Defaults

DEF_vfDomain = [0 2*pi];


% -- Check arguments

if (nargin < 2)
   help circ_ksdensityn;
   error('circ_ksdensityn:Usage', ...
         '*** circ_ksdensityn: Incorrect usage');
end

% - Do we have univariate data?
if (isvector(mfObservations))
   warning('circ_ksdensityn:univariate', ...
           '--- circ_ksdensityn: Assuming single dimensionality for vector data.');

	% - Reshape everything
   mfObservations = reshape(mfObservations, [], 1);
   mfPDFSamples = reshape(mfPDFSamples, [], 1);
end

vnObservationsSize = size(mfObservations);
vnPDFSamplesSize = size(mfPDFSamples);

if (~exist('mfDomains', 'var') || isempty(mfDomains))
   mfDomains = repmat(DEF_vfDomain, vnObservationsSize(2), 1);
end

% - Do we need to estimate fSigma?
if (~exist('vfSigmas', 'var'))
   % - Sigma will be estimated
   vfSigmas = [];
end

% - If weights are not provided, weight each observation equally
if (~exist('vfWeights', 'var'))
   vfWeights = 1;

% - Check the number of observations matches the number of weights
elseif (size(mfObservations, 1) ~= numel(vfWeights))
   error('circ_ksdensityn:Usage', ...
         '*** circ_ksdensityn: The number of observations must be equal to the number of weights.');
end

% - Check the dimensionality of all arguments
if (vnPDFSamplesSize(2) ~= vnObservationsSize(2)) || (size(mfDomains, 1) ~= vnObservationsSize(2)) 
   error('circ_ksdensityn:Usage', ...
         '*** circ_ksdensityn: The number of dimensions must match for ''msObservations'', ''mfPDFSamples'' and ''mfDomains''.');
end


% -- Clean data by removing NaNs

vbNanSample = any(isnan(mfObservations), 2);
mfObservations = mfObservations(~vbNanSample, :);
if (~isscalar(vfWeights))
   vfWeights = vfWeights(~vbNanSample);
end
vnObservationsSize = size(mfObservations);


% -- Map everything to [0..2pi]

vbPeriodicDomains = ~any(isnan(mfDomains) | isinf(mfDomains), 2);

for (nDim = find(vbPeriodicDomains)')
   % - Map observations to [0..2pi] and wrap
   mfObservations(:, nDim) = (mfObservations(:, nDim) - mfDomains(nDim, 1)) ./ diff(mfDomains(nDim, :)) .* 2*pi;
   mfObservations(:, nDim) = mod(mfObservations(:, nDim), 2*pi);
   
   % - Map PDF samples to [0..2pi] and wrap
   mfPDFSamples(:, nDim) = (mfPDFSamples(:, nDim) - mfDomains(nDim, 1)) ./ diff(mfDomains(nDim, :)) .* 2*pi;
   mfPDFSamples(:, nDim) = mod(mfPDFSamples(:, nDim), 2*pi);
end


% -- Wrap observations and PDF samples over periodic domains

mfRepObservations = mfObservations;
mfRepSamples = mfPDFSamples;
vbOrigSamples = true(size(mfRepSamples, 1), 1);
vfRepWeights = vfWeights;
for (nDim = find(vbPeriodicDomains)')
   % - Replicate observations over domain boundaries; replicate weights
   mfUpObservations = mfRepObservations;
   mfDownObservations = mfRepObservations;
   mfUpObservations(:, nDim) = mfRepObservations(:, nDim) + 2*pi;
   mfDownObservations(:, nDim) = mfRepObservations(:, nDim) - 2*pi;
   mfRepObservations = [mfRepObservations;
                        mfUpObservations;
                        mfDownObservations]; %#ok<AGROW>
	
	% - Replicate samples for estimating bin volumes
   mfUpSamples = mfRepSamples;
   mfDownSamples = mfRepSamples;
   mfUpSamples(:, nDim) = mfUpSamples(:, nDim) + 2*pi;
   mfDownSamples(:, nDim) = mfDownSamples(:, nDim) - 2*pi;
   mfRepSamples = [mfRepSamples;
                   mfUpSamples;
                   mfDownSamples]; %#ok<AGROW>
	vbOrigSamples = [vbOrigSamples;
                    false(size(mfRepSamples, 1)*2, 1)]; %#ok<AGROW>
                     
	% - Replicate weights, if required
   if (~isscalar(vfRepWeights))
      vfRepWeights = repmat(vfRepWeights, 3, 1);
   end
end


% -- Estimate sigmas, if necessary

if (isscalar(vfWeights))
   vfSTDWeights = ones(vnObservationsSize(1), 1);
else
   vfSTDWeights = vfWeights;
end

if (isempty(vfSigmas))
   parfor (nDim = 1:size(mfObservations, 2))
      if (vbPeriodicDomains(nDim))
         vfSigmas(nDim) = (4/3)^(1/5) * circ_std(mfObservations(:, nDim), vfSTDWeights) * (vnObservationsSize(1)^(-1/5)); %#ok<PFBNS>
      else
         vfSigmas(nDim) = (4 * nanstd(mfObservations(:, nDim), vfSTDWeights)^5 / 3 / vnObservationsSize(1))^(1/5);
      end
   end
end


% -- Perform kernel density estimate for each dimension in turn

parfor (nDim = 1:size(mfObservations, 2))
   if (isscalar(vfWeights))
      mfEstimate(:, nDim) = ksdensity(mfRepObservations(:, nDim), mfPDFSamples(:, nDim), 'width', vfSigmas(nDim));
   else
      mfEstimate(:, nDim) = ksdensity(mfRepObservations(:, nDim), mfPDFSamples(:, nDim), 'weights', vfRepWeights, 'width', vfSigmas(nDim));
   end
end

% - Form product over dimensions for independent multivariate estimate
vfEstimate = prod(mfEstimate, 2);


% -- Correct scaling of PDF estimate

% - Work out which version of triangulation classes to use
if (exist('delaunayTriangulation', 'class'))
   fhDelaunay = @delaunayTriangulation;
   fhVoronoi = @(dtRep, mfSamples)voronoiDiagram(dtRep);

elseif (exist('DelaunayTri', 'class'))
   fhDelaunay = @DelaunayTri;
   fhVoronoi = @(dtRep, mfSamples)voronoiDiagram(dtRep);
   
else
   % - Fall back to classic matlab functions
   fhDelaunay = @(x)[];
   fhVoronoi = @(dtRep, mfSamples)voronoin(mfSamples);
end
   
% - Compute a Voronoi partitioning of sample bins
dtTri = fhDelaunay(mfRepSamples);
[mfVoronoiVert, cvnVornoiCells] = fhVoronoi(dtTri, mfRepSamples);

% - Include Voronoi cells only for desired PDF samples
cvnVornoiCells = cvnVornoiCells(vbOrigSamples);

% - Estimate volumes for edge samples (i.e. cells connecting vertex 1 at infinity)
vnEdgeSamples = find(cellfun(@(c)(any(c == 1)), cvnVornoiCells));

% - Append edge samples as Voronoi vertices and replace infinity points for
% edge samples with original sample location
nStartAppendedVertices = size(mfVoronoiVert, 1);
mfVoronoiVert = [mfVoronoiVert; mfPDFSamples(vnEdgeSamples, :)];

for (nSample = 1:numel(vnEdgeSamples))
   % - Replace infinite vertex for this cell
   cvnVornoiCells{vnEdgeSamples(nSample)}(cvnVornoiCells{vnEdgeSamples(nSample)} == 1) = ...
      nSample + nStartAppendedVertices;
end

% - Estimate volumes for all cells
vfBinVol = cellfun(@(c)(bin_volume(mfVoronoiVert(c, :))), cvnVornoiCells);

% - Scale PDF estimate ideally
vfEstimate = vfEstimate ./ sum(vfEstimate .* vfBinVol) .* sum(vfWeights);

% --- END of circ_ksdensityn FUNCTION ---


% circ_std - FUNCTION Estimate the weighted circular standard deviation of a dataset
function [s] = circ_std(alpha, w)

% - Deal with scalar weight values
if (isscalar(w))
   w = repmat(w, size(alpha));
end

% - Compute mean resultant vector length
r = circ_r(alpha,w);

% - Estimte standard deviation
s = sqrt(2*(1-r));


% circ_r - FUNCTION Compute the weighted resultant of a dataset
function r = circ_r(alpha, w)

% - Compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),1);

% - Obtain vector length
r = abs(r)./sum(w,1);


% bin_volume - FUNCTION Return the volume of a bin using the convex hull
function fVol = bin_volume(mfBinVertices)

[~, fVol] = convhulln(mfBinVertices);

% --- END of circ_ksdensity.m ---


