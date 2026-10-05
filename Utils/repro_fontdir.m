function directory = repro_fontdir()
% Set SPECTRAL_FONT_DIR to a directory containing licensed Arial font files.
directory=getenv('SPECTRAL_FONT_DIR');
if isempty(directory) && ismac
    directory='/System/Library/Fonts/Supplemental';
elseif isempty(directory) && ispc
    directory=fullfile(getenv('WINDIR'),'Fonts');
end
assert(isfolder(directory), ...
    'Set SPECTRAL_FONT_DIR to the directory containing Arial.ttf and Arial Bold.ttf.');
end
