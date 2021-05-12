FROM opensciencegrid/osgvo-r:4.0.2

RUN mkdir /home/analysis

COPY installed_packages.rda /home/analysis 
COPY restore_packages.R /home/analysis

RUN Rscript /home/analysis/restore_packages.R
