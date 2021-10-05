# Template for scRNA visualizations with Shiny
![GitHub last commit](https://img.shields.io/github/last-commit/alemenze/magic_shiny_scRNA_template)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![made with Shiny](https://img.shields.io/badge/R-Shiny-blue)](https://shiny.rstudio.com/)

# Initial RData Processing
To host datasets, first generate an intial RData object containing the initial processing. This reduces the upfront computational burden to allow the focusing on visualization and exploration. 

## Running the App Locally
This Shiny App has been built in to a docker container for easy deployment. You can build the image yourself (and thereby customize any ports you need) after downloading it:
```
docker build -t scTemplate .
docker run -d --rm -p 8080:8080 scTemplate
```
And it should be hosted at localhost:8080

## Hosting datasets on GCloud
```
PROJECTID=$(gcloud config get-value project)
PROJECTNAME="ProjectName"
docker build . -t gcr.io/$PROJECTID/$PROJECTNAME
```
```
docker push gcr.io/$PROJECTID/$PROJECTNAME
```

See Deployment. 