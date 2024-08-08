# CosmologyMaskAnalyser

This is a library designed to analyse cosmological masks and catalogues. The library itself is in the folder named "toolkit" - this needs to be installed in the venv environment in order to run the code (in the "code" folder) or the batch files (in the "batch_files") folder.

## Features

- toolkit.toolkit:
-- This adds classes for both Healpix and Pixell masks, providing file handling, plotting scripts and vectorised querying for both.
-- This adds a class for cluster catalogues, providing file handling and plotting.
-- The bin map classes implement various tesselation algorithms (primarily Healpix).
- toollit.weights:
-- This provides various weighting algorithms used for the calculations.
- toolkit.data:
-- This is a wrapper for various data files - the currently referenced files are not included on the repository, this file needs to be updated for the code to be used.
