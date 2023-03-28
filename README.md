# Fresh-API 2: Development of a data science machine learning platform for food quality
**Background** : Approximately 88 million tonnes of food are wasted annually across the European Union due to spoilage, incurring serious financial and environmental implications. The safety of poultry-containing food products can be compromised when the meat is prepared and packaged. The presence of pseudomonads in aerobic conditions can produce unpleasant slime and off-odours compromising the quality of the meat. European Commission guidelines dictate that the quality of fresh meat be determined by total viable counts (TVC) of bacteria or counts of Enterobacteriaceae. Obtaining this type of data from standard colony-counting methods for meat samples is time-consuming and labor intensive with retrospective results. Multispectral imaging (MSI), Fourier Transform Infrared (FTIR) spectroscopy, Electronic nose (eNose), High Performance Liquid Chromatography (HPLC) and Gas Chromatography Mass Spectrometry (GCMS) technologies are being widely used to provide an alternative way to assess food quality as rapid and limited or non-invasive methods

**Aims and objectives** : To create a scalable platform based on Representational State Transfer (REST) and application of APIs (Application Programming Interface) for the rapid quality assessment of chicken-thigh fillet, chicken-burger samples and chicken-breast fillet. This will be achieved through the implementation of machine learning regression algorithms based on MSI, FTIR, HPLC, GCMS, eNose data to predict microbial loads and classification algorithms to predict sensory scores.

**Approach** : The approach includes the employment of multiple machine learning algorithms to achieve the optimum prediction accuracies. In particular, five algorithms have been implemented: Linear Regression, K-Nearest Neighbor, Random Forest, Support VectorMachines with Polynomial Kernel, Radial Basis Function Kernel and naïve bayes. Analyses and modelling were conducted via R language and environment for statistical computing and graphics, including the API (Plumber) which acts as a “pipeline” facilitating data requests and returning a response as a predicted value of the bacterial count or sensory score of a given sample 

## Table of Contents
* [Fresh-API](#Fresh-API)
* [Installation](#Installation)
  * [Dependencies](#Dependencies)
  * [Install Fresh-API from source](#Install-Fresh-API-from-source)
  * [Step 1 - Installing Docker](#Step-1---Installing-Docker)
  * [Step 2 - Building the Docker Image](#Step-2---Building-the-Docker-Image)
  * [Step 3 - Running the Docker container](#Step-3---Running-the-Docker-container)
  * [Step 4 - Sending request to REST-API](#Step-4---Sending-request-to-REST-API)
    * [General Argument Description](#General-Argument-Description)
    * [Sending request to prediction and classification endpoints](#Sending-request-to-prediction-and-classification-endpoints)
    * [Sending request to data inspection and visualization endpoints](#Sending-request-to-data-inspection-and-visualization-endpoints)
    * [Sending request to model training endpoints](#Sending-request-to-model-training-endpoints)
* [JSON file format specification](#JSON-file-format-specification)
* [Additional commands for Docker](#Additional-commands-for-Docker)
* [Troubleshooting](#Troubleshooting)
  * [HTTP General Error codes](#HTTP-General-Error-codes)
  * [Possible type of errors in model serving endpoints](#Possible-type-of-errors-in-model-serving-endpoints)
  * [Possible errors in data preparation and model training endpoints](#Possible-errors-in-data-preparation-and-model-training-endpoints)

## Fresh-API 
The purpose of the application is to provide a platform to predict the total bacterial counts present on a given meat sample through the implementation of REST API. Fresh-API is build upon the technology of RESTful web based API consisting of an **Endpoint URL** and **Body Data** which can be transmitted via **HTTP request methods** by sending a single JSON-encoded data string. Machine Learning regression models are developed onto the platform and validated against new data coming to the system by sending the prediction results through the REST backend.

# Installation
## **Dependencies**
 * UNIX Operating Systems  
 * [Docker](https://www.docker.com/why-docker)
 * R (tested on version 3.6.3)
 * The following R libraries:
 
      | Package | Version | Package | Version |
      |-------|-------|--------|--------|
      | `plumber` | 0.4.6 | `ggplot2` | 3.3.0|   
      | `rmarkdown` | 2.1 | `ggpubr` | 0.3.0 |   
      | `caret` | 6.0.86 | `zoo` | 1.8.7|
      | `randomForest` | 4.6.14 | `factoextra` | 1.0.7| 
      | `tinytex` | 0.21 | `gmodels` | 2.18.1 |
      | `jsonlite` | 1.6.1 | `dplyr` | 0.8.5|   
      | `openxlsx` | 4.1.4 | `assertthat` | 0.2.1|
      | `stats` | 4.1.4 |  | |      
      | `prospectr` | 0.2.0 | | |
      | `kernlab` | 0.9.29 | | |   
      | `readxl` | 1.3.1 | | |
 
## **Install Fresh-API from source**  
```
$ git clone https://github.com/FadyMohareb/freshAPI_GP.git

$ cd freshAPI_GP/
```

## **Step 1 - Installing Docker**

To download the latest version, install Docker from the official Docker repository. This section guides you how to do that.

First, in order to ensure the downloads are valid, add the GPG key for the official Docker repository to your system:
```
$ curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
```
Add the Docker repository to APT sources:
```
$ sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
```
Update the package database with the Docker packages from the newly added repository:
```
$ sudo apt-get update
```
To ensures that you are running the installation from the official Docker repository, run the command:
```
$ apt-cache policy docker-ce
```
This should give an output similar to the following:
```
docker-ce:
  Installed: 5:19.03.8~3-0~ubuntu-xenial
  Candidate: 5:19.03.8~3-0~ubuntu-xenial
  Version table:
 *** 5:19.03.8~3-0~ubuntu-xenial 500
        500 https://download.docker.com/linux/ubuntu xenial/stable amd64 Packages
        100 /var/lib/dpkg/status
     5:19.03.7~3-0~ubuntu-xenial 500
        500 https://download.docker.com/linux/ubuntu xenial/stable amd64 Packages
 ```
From the output, you will notice that the docker-ce is not yet installed. However, the output will show the target operating system and the version number of the Docker. Please note that version numbers may differ depending on the time of installation.

Use the following command to install Docker:
```
sudo apt install docker-ce
```
This will install Docker, start the daemon and enable it to automatically start on boot. To confirm that the Docker is active and working, run:
```
sudo systemctl status docker
```
The command will provide the following output:
```
● docker.service - Docker Application Container Engine
   Loaded: loaded (/lib/systemd/system/docker.service; enabled; vendor preset: enabled)
   Active: active (running) since Mon 2020-03-23 12:34:27 GMT; 4 weeks 1 days ago
     Docs: https://docs.docker.com
 Main PID: 11575 (dockerd)
    Tasks: 31
   Memory: 1.9G
      CPU: 33min 41.448s
   CGroup: /system.slice/docker.service
           └─11575 /usr/bin/dockerd -H fd:// --containerd=/run/containerd/containerd.sock
```           
## **Step 2 - Building the Docker Image**
`docker` consists of chain of options and commands followed by arguments. The syntax takes this form:
```
$ docker [option] [command] [arguments]
```
To view all available subcommands, type:
```
$ docker
```
The `docker build` command builds Docker images from a Dockerfile and a “context”. A build’s context is the set of files located in the specified `PATH`:
```
$ docker build -t freshapi .
```
NOTE: `-t` or `--tag` is an option used for label a docker image in the `name:tag` format. In this case, name will be `freshapi`. The `.` operator specifies the current directory containing the 'Dockerfile'. Optionally, the version of an image can be tagged in which the user would like to run the container. For example, `docker build -t freshapi:2.0 .` 

## **Step 3 - Running the Docker container**
The `docker run` command creates a container from a given image and starts the container using a given command or entrypoint. To run the resulting docker image, the command is as follows:
```
$ docker run -v $(pwd):/freshAPI_GP --name [container name of choice] -d freshapi
```
NOTE: The `-v` flag mounts the current working directory into the container. `pwd` is a command that "print current/working directory". Putting the command in `$()` takes the command within parenthesis, runs it in a subshell and yields back the absolute path to our project directory.`-d` or `--detach` is an option for running the container in background. `--name` is used to assign a name to the container. 

In order to execute commands on running containers, `docker exec` is used to specify the container name as well as the command to be executed on this container:
```
$ docker exec -it [container name of choice] bash
```
NOTE: The combination of `-i` and `-t` allows the container to run in an interative mode and provides access to the terminal to execute further commands in the application. The `bash` command creates a new Bash session for the container.

## **Step 4 - Sending request to REST-API**
The design and features of Fresh-API is to accept request containing data and query parameters defined by the developers. curl is a command line tool for transferring and receiving HTTP requests and responses. The description of all endpoints argument used are as follows:
## **General Argument Description**
| Argument | Description |
|----------|-------------|
| `filename` | Is the file name with its extension, which contains data to be processed or used, in this case it can be `.csv`, `.xlsx` or `.rds` file. |
| `filename1` | 	Is the file name with its extension, which contains data to be processed or used, in this case it can be `.csv`, `.xlsx` file only. |
| `filename2` | Is the file name with its extension, which contains data to be processed or used, in this case it can be `.csv`, `.xlsx` file only. |
| `path`	| Is a folder name in which the data are going to be read from, or path after the `root` working directory. | 
| `platform` | Is the name of analytical platform that has been used to generate the data, i.e `HPLC`, `GCMS`, `eNose`, `MSI` or `FTIR` only. |
| `product`	| Is the name of meat product type, of which the data has been obtained. Here, in this case chicken burger, chicken thigh fillet or chicken breast fillet denoted by `CBG`, `CTF` or `CBF` respectively. | 
| `products` | Is the name of meat product type, of which the data has been obtained. Here, in this case chicken burger, chicken thigh fillet or chicken breast fillet denoted by Burger, Thigh or Breast respectively. |
| `bacteria`	| Is the name of bacterial type, i.e `TVC`, `Ps`, `Bth` or `LAB`. |
| `replicates`	| Is the number of replicate samples present in the data, it’s a positive numeric values starting from `1` |
| `Rm`	| It is a logical argument used to remove or impute missing values in the dataset, it’s always yes to remove or no to impute. |
| `PC`	| Is used to choose which principal component to plot i.e `PC12` denote plotting Pc1 vs Pc2 or `PC13` plots Pc1 vs Pc3. Only this `PC12` or `PC13` can be used. |
| `iteration`	| Is a positive numeric value passed to specify number of iteration to train over, that start from `1`|
| `Batches` | If data for machine learning contains batches and user want to train on some batches and test on other. To train on mixed batch use `mixed` or to train on batch use for example `R1;R2`, which is `R1` for training and `R2` for testing. |
| `pram`	| This argument is used to specify [`hyper-parameters`](https://topepo.github.io/caret/random-hyperparameter-search.html) for [`tune grids`](https://rpubs.com/Mentors_Ubiqum/tunegrid_tunelength) and [`grid search`](https://topepo.github.io/caret/random-hyperparameter-search.html), depending on the ML algorithm used this argument can accommodate different [tuning](https://topepo.github.io/caret/train-models-by-tag.html)(caret package) parameters on it. Example for `svmPoly` tuning parameters can be:  `pram=C~1,0.2,0.3,0.01,1;degree~1,2,3;scale~0.1,0.5`. The same way tuning parameters can be passed to other ML algorithms. Note: the tilde Symbol `~` is used to assign values to parameters and semicolon symbol `;` is used to separate parameters. |
| `method` | Passes machine learning algorithm to be used to train(`rf`, `svmLinear`, `knn`, `svmPoly`, `svmRadial`) |
| `round` | is a logical argument which accepts `yes` or `no`, `yes` will be incase of rounding FTIR wavelength data, and no is used when .rds object needs to be converted to .csv |
| [`tunelength`](https://rpubs.com/Mentors_Ubiqum/tunegrid_tunelength) | Passes numeric value that start from 0, Is used to specify to how many times to try different default values for the main parameter. 
## **Sending request to prediction and classification endpoints**
#### **Bacterial count prediction**
##### Eg.
```
curl --data data.json "localhost:8080/Predict?platform=eNose&product=CBF&bacteria=TVC"
```
#### **Sensory score prediction**
##### Eg.
```
curl --data data.json "localhost:8080/classify?platform=eNose&product=CBF"
```
## **Sending request to data inspection and visualization endpoints**
#### **Inspecting missing values in dataset**
##### Eg.
```
curl --data null "localhost:8080/row_index_NA?filename=enose_breastAUA.xlsx&path=Demo"
```
#### **Dataset missing value handling**
##### Eg.
```
curl --output enoseN.csv "localhost:8080/rm-im-NA?platform=eNose&path=Demo&filename=enose_breastAUA.csv&Rm=No&replicates=1"
```
#### **Replicate averaging in dataset**
##### Eg.
```
curl --output eNose_rep.csv "localhost:8080/replicates_tre?platform=eNose&replicates=3&filename=enoseN.csv&path=De&products=CBF"
```
#### **PCA plot generation based on batch sample present in dataset**
##### Eg.
```
curl --output FL_PC12.png "localhost:8080/PC?platform=FTIR&filename=chickenBurger_FTIR_raw_woo_v2.xlsx&path=Demo&PC=PC12”
```
#### **Replicate arrangement and quality checking**
##### Eg.
```
curl --data null "localhost:8080/check_rep?filename=enose.csv&path=Demo&replicates=3”
```
#### **Combining new datasets and existing dataset for machine learning**
##### Eg.
```
curl --data newd.csv "localhost:8080/Combine_dataset?platform=FTIR&filename1=FTIR1.csv&filename2=FTIR2.csv&path=Demo”
```
## **Sending request to model training endpoints**
#### **Training machine learning model**
##### Eg.
```
curl "localhost:8080/MLearning?platform=HPLC&bacteria=TVC&Batches=mixed&method=svmLinear&iteration=10&path=De&filename=HPLC_rep.csv&pram=C~0.1,0.2,0.3,0.4,1&tunelength=0"
```
#### **Plotting model summary plots from rds object to PNG**
##### Eg.
```
curl --data AP.png "localhost:8080/plots?filename=sp.rds&path=Demo&model_plot=Pre_Act”
```
#### **Generating CSV file of predicted and observed bacterial count from RDS object**
##### Eg.
```
curl --data AP.csv "localhost:8080/Predicted_Actual?filename=dm.rds&path=Demo&platform=HPLC&round=no”
```
#### **Rounding wavelength of FTIR data**
##### Eg.
```
curl --data FTIR.csv "localhost:8080/Predicted_Actual?filename=FTIRd.csv&path=Demo&platform=FTIR&round=yes”
```
The following table provides a list of endpoints and their corresponding query parameters:
   |Endpoints|Query parameters| 
   |---------|----------------|
   |`#* @post /Predict`|`product= CTF, CBG, CBF`; `bacteria= TVC, Ps, Bth, LAB`; `platform=MSI, FTIR, GCMS, HPLC, eNose`|
   |`#* @post /classify`|`product= CBF` ; `bacteria= TVC` |
   |`#* @post /row_index_NA`| `filename=..` ; `path=...` |
   |`#* @get /rm-im-NA`|`platform= MSI, FTIR, eNose, GCMS, HPLC` ; `filename= ...`; `path= ...` ; `replicates= 1...` ; `Rm= yes, no`|
   |`#* @get /replicates_tre`|`platform= MSI, FTIR, eNose, GCMS, HPLC` ; `product= CBG, CTF, CBF` ; `replicates= 1...` ; `filename=...` ; `path=...` |
   |`#* @get /PC`|`platform= MSI, FTIR, eNose, GCMS, HPLC` ; `PC= PC12, PC13` ; `path=...` ; `filename= ...`|
   |`#* @post /check_rep`|`filename=...` ; `path=...` ; `replicates= 1...`|
   |`#* @get /Combine_dataset`|`platform= MSI, FTIR, eNose, GCMS, HPLC` ; `filename1=...` ; `filename1=...` ; `path=...`|
   |`#* @get /MLearning`|`platform= MSI, FTIR, eNose, GCMS, HPLC` ; `bacteria= TVC, Ps, Bth, LAB`; `Batches=mixed, R1;R2, ...` ; `method=svmLinear, svmPoly, svmRadial, knn, rf` ; `iteration=1...` ; `filename=...` ; `path=...` ; `pram= ...` ; `tunelength=0, 1, ...`|
   |`#* @get /ML_plots`|`model_plot= Iter, Pre_Act` ; `filename=...` ; `path=...` |
   |`#* @get /Predicted_Actual`|`platform= MSI, FTIR, eNose, GCMS, HPLC` ; `filename=...` ; `path=...` ; `round=yes, no`|
   

### Legend
1)	`platform` column contains information about the type of platform:
    *	`FTIR` = Fourier Transformed Infrared spectroscopy
    *	`MSI` = Multi-spectral iamging.
    * `HPLC` = High Performance Liquid Chromatography
    * `GCMS` = Gass Chromatography Mass Spectrometry
    * `eNose` = Electronic nose
2)	`product` = type of product:
    *	`CB` = chicken burger
    *	`CTF` = chicken thigh fillet
    * `CBF` = chicken breast fillet
4)	`bacteria` = type of bacterial counts:
    *	`TVC` = total viable counts
    *	`Ps` = Pseudomonas spp.
    *	`Bth` = Brochothrix thermosphacta
    *	`LAB` = lactic acid bacteria
5)	`model` = machine learning method based on what the model was built:
    *	`knn` = k-nearest neighbors algorithm (k-NN)
    *	`rf` = random forest algorithm
    *	`svmLinear` = support-vector machine (SVM) with linear kernel
    *	`svmRadial` = support-vector machine (SVM) with radial kernel
    *	`svmPoly` = support-vector machine (SVM) with polynomial kernel
    * `lm` = linear regression (linear model)

NOTE: `--data` or `-d` denotes the `curl` command used for passing data to the request body and `--output` or `-o` denotes the file name of output object. `platform`, `product`, `bacteria` and `model` parameters were passed to the endpoint using “query strings”. The (`?`) appended to the URL indicates the start of the query string. In the query string, each parameter is concatenated with other parameters through the ampersand (`&`) symbol.

# JSON file format specification
* In the case of predicting bacterial counts, data derived from an analytical platform such as MSI should contain the 18 mean values as features at the beginning of a JSON file.
* For FTIR derived samples the file should contain wavelengths in the range of 1001-4000 nm. 
* For eNose acquired data it should containe 12 features with exact sensory array names (see eNose data).
* For HPLC acquired data it should containe 9 features with exact biochemical names (see HPLC data)
* For GCMS acquired data it should containe 41 features with exact biochemical names or volatile metabolites names (see GCMS data)

# Additional commands for Docker
   | USAGE | Commands |
   |-------|----------|
   |Display logs of a container|`$ docker logs [container name]`|
   |List all existing containers (running and not running)|`$ docker ps -a`|
   |List your images|`$ docker image ls`|
   |Stop a specific container|`$ docker stop [container name]`|
   |Stop all running containers|`$ docker stop $(docker ps -a -q)`|
   |Delete a specific container (only if stopped)|`$ docker rm [container name]`|
   |Delete all containers (only if stopped)|`$ docker rm $(docker ps -a -q)`|
   |Delete all unused containers, unused images and networks|`$ docker system prune -a --volumes`|
# Troubleshooting
###### **HTTP General Error codes**
| Http status code  | Explanation |
|-------|----------|
|  200 Success     | Request is successful |
|  400 Bad Request | HTTP request sent to the server has invalid syntax          |
|  404 Not Found    | unable to locate the requested file or resource, mainly these is encounterd when the url is incorrect or file specified doesn't exist |
|  500 Internal Server Error   | Request cannot be processed, these may happen when argument passed are not proper (see docker logs for detail)|



###### **Possible type of errors in model serving endpoints**

```
Error: No data provided
```
```
Error Occured: Please choose correct platform: FTIR | HPLC | eNose | GCMS | MSI 
```
```
Error occured: Please make sure product type is correct (CBG or CFT or CBF)
```
```
Error occured: Please make sure bacteria type is correct (TVC or Ps or LAB or Bth )
```
###### **Possible errors in data preparation and model training endpoints**
Note: Here, most the endpoint are serializer content-type. Therefore, it's not possible to get error message sometimes, however check output file if it contains the desired output. And use the command below to check the `logs` details in the `docker container` which gives the below information about the error.
```
$ docker logs [container name]
```
```
Error occured: Please check your file it should contain .xlsx or .csv
```
```
Error Occured: Please check your file contains appropraite replicate number
```
```
Error Occured: please make sure you choose correct platform: FTIR | MSI | HPLC  | GCMS  | eNose
```
```
Error occured: Please check PC argument has to be PC13 or PC12
```
```
Error occured: Please make sure tunelength is zero or greater
```
```
Error occured: Please make sure iteration is 1 to 100 values only
```
```
Error occured: Please choose corret machine learning algorithm: svmRadial | svmPoly | svmLinear | knn | rf
```
```
Error occured: Please choose correct model_plot argument: Iter or Pre_Act 
```
```
Error occured: Please choose correct filename: sp.rds or st.rds
```
```
Error Occured: please choose correct file name with .rds extention
```
```
Warning message: Please check some of your samples doesn't contain replicate of :
```
```
Warning message: Please check your samples replicates aren't underneath of each other 1 : n or some of your samples doesn't contain replicate of n
```
```
Warning message: unable to check replicates
```
