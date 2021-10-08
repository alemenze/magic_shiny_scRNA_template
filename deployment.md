# Deployment notes

```
PROJECTID=$(gcloud config get-value project)
PROJECTNAME='projectname'
cluster_name=
image_name=
service_name=

docker build . -t gcr.io/$PROJECTID/$PROJECTNAME
```
```
docker push gcr.io/$PROJECTID/$PROJECTNAME
```

Since these will exceed what we can use within cloudrun, it would require a kubernetes cluster or similar:
## Starting the GKE cluster
```
gcloud container clusters create $cluster_name --machine-type=e2-standard-4 --enable-autoscaling --min-nodes=1  --max-nodes=3 --num-nodes=1 --zone=us-east4-a

```
*Adjust the machine as needed to fit the size of project*

Get your credentials to the cluster
```
gcloud container clusters get-credentials $cluster_name --zone us-east4-a
```

## Adding a project to the cluster
Start here if you already have a live cluster. Deploy the image to the cluster
```
kubectl create deployment $image_name --image=gcr.io/$PROJECTID/$PROJECTNAME
```

Set up the scaling. You can tweak this to specific projects as needed.
```
kubectl scale deployment $image_name --replicas=1
kubectl autoscale deployment $image_name --cpu-percent=80 --min=1 --max=5
```

You may have to wait a few minutes for GKE to build your pod- you can monitor it with
```
kubectl get pods
```

Once you see the pod running, you need to expose the ports. Change this to fit the specific in ports and target ports
```
kubectl expose deployment $image_name --name=$service_name --type=LoadBalancer --port 80 --target-port 8080
```
*might need to update this to fit a nginx controller*

Get the actual external IP address
```
kubectl get service
```