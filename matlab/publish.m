
% Options
publishOptions.format = 'html';
publishOptions.showCode = false;
publishOptions.evalCode = false;

% Generate
publish('README.m', publishOptions);

% Show
web('html/README.html')