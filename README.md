# preseimei
Pre-Processing for raw fits. Dark & flat subtraction. Give WCS information using astrometry.net. 

## Requirement 
You need to install astrometry.net beforehand http://astrometry.net/doc/readme.html  
Confirm that you can use the command "solve-field". 

Then, you need to get *index files* from http://data.astrometry.net/4200/, and put them to data directory for astrometry.net (possibly /usr/local/astrometry/data)

Read http://data.astrometry.net/4200/ carefully. Index files are needed for astrometry.net to solve WCS. 


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
