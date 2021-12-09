#!/bin/bash

echo "--> Initializing project to push as shiny app"

# Get the path and set up config
echo "--> Reading in your config..."
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

configFile="${DIR}/config.json"
CLUSTER_NAME=`jq -r '.kubernetes .cluster_name' ${configFile}`
MACHINE_TYPE=`jq -r '.kubernetes .machine_type' ${configFile}`
NODE_COUNT=`jq -r '.kubernetes .node_count' ${configFile}`
ZONE=`jq -r '.kubernetes .zone' ${configFile}`
PROJECTID=`jq -r '.kubernetes .PROJECTID' ${configFile}`
PROJECTNAME=`jq -r '.shiny .PROJECTNAME' ${configFile}`
IMAGE_NAME=`jq -r '.shiny .image_name' ${configFile}`
SERVICE_NAME=`jq -r '.shiny .service_name' ${configFile}`

# Check if the cluster is operational
echo "--> Checking on the cluster..."
active_clusters="$(gcloud container clusters list | awk '{print $1}' )"
cluster_alive=false

for item in $active_clusters
do
    if [ "$CLUSTER_NAME" == "$item" ]; then
        cluster_alive=true
        echo "--> Cluster $CLUSTER_NAME is currently active, preparing to deploy image..."
    fi
done

# Build cluster if not
if [ "$cluster_alive" = false ]; then
    echo '--> Preparing to build the cluster..'
    gcloud container clusters create ${CLUSTER_NAME} \
        --machine-type ${MACHINE_TYPE} \
        --enable-autoscaling \
        --num-nodes ${NODE_COUNT} \
        --min-nodes=1 \
        --max-nodes=10 \
        --zone ${ZONE} \

    nodecount="$(kubectl get node | awk '{print $2}' | grep Ready | wc -l)"
    while [[ ${nodecount} -lt ${NODE_COUNT} ]] ; do echo -n $(date) ; echo " : ${nodecount} of ${NODE_COUNT} nodes ready" ; sleep 15 ; nodecount="$(kubectl get node | awk '{print $2}' | grep Ready | wc -l)" ; done
    echo "--> Cluster $CLUSTER_NAME has been built!"
fi

# Build your image
echo "--> Building your docker image..."
docker build . -t gcr.io/$PROJECTID/$PROJECTNAME --quiet
echo "--> Image built!"

# Deploy
echo "--> Preparing to deploy the image..."
docker push gcr.io/$PROJECTID/$PROJECTNAME
echo "--> Image pushed!"
gcloud container clusters get-credentials ${CLUSTER_NAME} --zone us-east4-a

#Check for existing image
echo "--> Checking for active deployment..."
deployment_alive=false
active_deploy="$(kubectl get deployment $IMAGE_NAME | awk '{print $1}')"
for item in $active_deploy
do
    if [ "$IMAGE_NAME" == "$item" ]; then
        deployment_alive=true
        echo "--> $IMAGE_NAME is currently active, preparing to restart with new image..."
    fi
done

if [ "$deployment_alive" = false ]; then
    echo "--> Preparing to start the deployment..."
    kubectl create deployment $IMAGE_NAME --image=gcr.io/$PROJECTID/$PROJECTNAME
    kubectl scale deployment $IMAGE_NAME --replicas=1
    kubectl autoscale deployment $IMAGE_NAME --cpu-percent=80 --min=1 --max=5

    podStatus="$(kubectl get pods --selector=app=${IMAGE_NAME} | grep ^${IMAGE_NAME} | awk '{print $3}')"
    while [[ ! x${podStatus} == xRunning ]] ; do echo -n $(date) ; echo " : pod status : ${podStatus} " ; sleep 30 ; podStatus="$(kubectl get pods --selector=app=${IMAGE_NAME} | grep ^${IMAGE_NAME} | awk '{print $3}')" ; done

    kubectl expose deployment $IMAGE_NAME --name=$SERVICE_NAME --type=LoadBalancer --port 80 --target-port 8080
    echo "--> Deployment ready!"

else
    echo "--> Preparing to restart the deployment..."
    kubectl rollout restart deployment $IMAGE_NAME
    podStatus="$(kubectl get pods --selector=app=${IMAGE_NAME} | grep ^${IMAGE_NAME} | awk '{print $3}')"
    while [[ ! x${podStatus} == xRunning ]] ; do echo -n $(date) ; echo " : pod status : ${podStatus} " ; sleep 30 ; podStatus="$(kubectl get pods --selector=app=${IMAGE_NAME} | grep ^${IMAGE_NAME} | awk '{print $3}')" ; done
    echo "--> Deployment ready!"
    
fi

# Print out the service IP address
SERVICE_IP=`kubectl get svc ${SERVICE_NAME} | awk '{ print $4}' | tail -n 1`
echo $"${SERVICE_NAME} IP: ${SERVICE_IP}"
while [ "${SERVICE_IP}" = '<pending>' ] || [ "${SERVICE_IP}" = "" ]
do
    echo "Sleeping 30s before checking again"
    sleep 30
    SERVICE_IP=`kubectl get svc ${SERVICE_NAME} | awk '{ print $4}' | tail -n 1`
    echo $"${SERVICE_NAME} IP: ${SERVICE_IP}"
done