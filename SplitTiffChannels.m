function SplitTiffChannels

[ filename, pathname ] = uigetfile( '*.tif', 'Pick multichannel tiff' );
fullpath=[pathname filename];
%tmpInfo = imfinfo(fullpath);

[header,Aout,cmap] = scim_openTif(fullpath, 'channels', 2, 'write', [ strrep( fullpath, '.tif', '_channel2' )]);
[header,Aout,cmap] = scim_openTif(fullpath, 'channels', 3, 'write', [ strrep( fullpath, '.tif', '_channel3' )]);

end