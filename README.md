# diHMM v1.0
Here presents an improved implementation of diHMM, contributed by Stephanos Tsoucas and Yan Kai.

diHMM is a novel chromatin segmentation for annotating chromatin states at two length scales. It was originally developed and described in the publication *Multi-scale chromatin state annotation using a hierarchical hidden Markov model* ([Marco et.al 2017](https://www.nature.com/articles/ncomms15011)). Detailed information of the project can be found at [this page](https://github.com/gcyuan/diHMM).

We applied the improved diHMM algorithm here to generate the multi-scale chromatin state annotations for the 127 human reference epigenomes in the [Roadmap and ENCODE consortia](http://www.roadmapepigenomics.org). Detailed information of the information on the 127 epigenomes can be found at [this site](https://egg2.wustl.edu/roadmap/web_portal/meta.html).

## Accessing the multi-scale chromatin state maps in the 127 epigenomes
We generated the chromatin state maps at the nucleosome (200bp resolution) and domain (4kb resolution) level. These maps are housed at *https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/*, and they can be download in the following way in terminal:
```
# specify the ID of the epigenome that you're interested in
EID=E001
# specify which level of state map you'd like to download. Could be "nucleosome" or "domain"
Type=domain
wget https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/$EID_$Type.bed.gz
```
The information on epigenome ID can be found [here](https://egg2.wustl.edu/roadmap/web_portal/meta.html).

## Installation
Go into the build dir and run
```
cmake ..
make
```
Then you can open a Python shell in the same dir and do
```
>>> import greet_ext
>>> greet_ext.greet()
'Hello world'
```

## Training a model
Training a diHMM model can be done by using the script *train.py*, after making necessary changes to input data path or other parameters.

## Applying a diHMM model for chromatin state annotation
Annotation can be done with the script *annotation.py*, after making necessary changes to input data path or other parameters.
