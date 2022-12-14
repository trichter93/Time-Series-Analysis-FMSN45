%
% Data file: dataCaffeine.m
%
% The file contains a time series relating to caffeine levels in instant 
% coffee. The data come from a cyclic process with period five units.
%
% This data was analyzed in 
%   Hamilton and Watt, "Interpreting partial autocorrelation functions of 
%   seasonal time series model", Biometrika 65: 135-140, 1978
% 
%   Chatfield, "Inverse autocorrelation", JRSS A 142:363-377, 1979.
%
%   Tiao and Tsay, "Identification of nonstationary and stationary ARMA
%   models", 1981 ASA Proceedings, 308-312.
%

% Reference: 
%   "An Introduction to Time Series Modeling", 4th ed, by Andreas Jakobsson
%   Studentlitteratur, 2021
%
data = [ 0.429 0.443 0.451 0.455 0.440 0.433 0.423 0.412 0.411 0.426 0.436 0.441 0.446 0.443 0.437 0.426 0.427 0.431 0.436 0.432  0.433 0.424 0.420 0.416 0.405 0.408 0.414 0.417 0.400 0.404 0.396 0.393 0.389 0.409 0.409 0.413 0.403 0.402 0.389 0.392  0.383 0.389 0.386 0.395 0.400 0.412 0.411 0.414 0.412 0.404 0.401 0.406 0.407 0.410 0.410 0.409 0.407 0.407 0.402 0.409  0.403 0.398 0.397 0.396 0.392 0.391 0.395 0.391 0.387 0.389 0.399 0.396 0.392 0.391 0.388 0.395 0.405 0.414 0.425 0.433  0.423 0.417 0.432 0.430 0.418 0.420 0.418 0.404 0.403 0.418 0.419 0.417 0.418 0.420 0.417 0.421 0.422 0.428 0.429 0.428  0.421 0.422 0.412 0.405 0.405 0.411 0.407 0.411 0.413 0.409 0.405 0.408 0.402 0.398 0.396 0.392 0.391 0.399 0.407 0.406  0.393 0.387 0.386 0.383 0.389 0.393 0.396 0.387 0.380 0.359 0.361 0.375 0.399 0.406 0.417 0.421 0.407 0.384 0.393 0.410  0.409 0.413 0.410 0.398 0.382 0.381 0.376 0.389 0.395 0.397 0.401 0.404 0.401 0.393 0.401 0.399 0.400 0.402 0.399 0.397   0.398 0.399 0.388 0.389 0.379 0.384 0.394 0.402 0.395 0.398 0.391 0.375 0.375 0.393 0.394 0.389 0.391 0.385 ];