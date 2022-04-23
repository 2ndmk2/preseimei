# preseimei
Pre-Processing for raw fits (SEIMEI TriCCS data in mind). Dark & flat subtraction. Give WCS information using astrometry.net. 

## Requirement 
You need to install astrometry.net beforehand http://astrometry.net/doc/readme.html  
Confirm that you can use the command "solve-field". 

## Get index files for astrometry.net
**After installing, you have to read http://astrometry.net/doc/readme.html carefully.**  

You need to get *index files* from http://data.astrometry.net/4200/, and  
put them to data directory for astrometry.net (possibly /usr/local/astrometry/data)  
Be careful that the total size can be ~30GB. 

Necessary index files depend on FOV, and it would be OK to get only index-4201~4205 files for TriCCS (6 * 12 arcmins^2).  
You can use the function for downloading files:  

from preseimei.install import download_index_files

## Install 
For normal install, 
* python setup.py install
* *For conda/pip env*, go to parental directory of tesarot, and "pip install ./preseimei"

## Install for development
For developers, 

* python setup.py develop
*  *For conda/pip env*, go to parental directory of tesarot, and "pip install -e ./preseimei"

## Test run
Go to tests folder
