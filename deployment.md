# Deployment notes
Based on [rand3k](https://github.com/randy3k/shiny-cloudrun-demo) repo

```
PROJECTID=$(gcloud config get-value project)
PROJECTNAME='projectname'
docker build . -t gcr.io/$PROJECTID/$PROJECTNAME
```
```
docker push gcr.io/$PROJECTID/$PROJECTNAME
```
```
gcloud run deploy --image gcr.io/$PROJECTID/$PROJECTNAME --platform managed --max-instances 1
```
Manually adjust CPUs and RAM applied to the container as it may be custom. 