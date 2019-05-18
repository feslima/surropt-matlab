% Get samples inside a hypercube bound from a larger sample set. Necessary
% for sampling set of Caballero methodology
%
% Inputs:
%       - xlhs: Original sample set set
%       - hlb:  HyperCube lower bound
%       - hub:  HyperCube upper bound
function [samples,fobssampled,gobsSampled] = getSamples(xlhs,fobs,gobs,hlb,hub)

xIndex = find(all(bsxfun(@lt,xlhs,hub) & bsxfun(@gt,xlhs,hlb),2));

samples = xlhs(xIndex,:);
fobssampled = fobs(xIndex);
gobsSampled = gobs(xIndex,:);