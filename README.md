# diHMM v1.0
Here presents an improved implementation of diHMM, contributed by Stephanos Tsoucas and Yan Kai.

diHMM is a novel chromatin segmentation for annotating chromatin states at two length scales. It was originally developed and described in the publication *Multi-scale chromatin state annotation using a hierarchical hidden Markov model* ([Marco et.al 2017](https://www.nature.com/articles/ncomms15011)). Detailed information of the project can be found at [this page](https://github.com/gcyuan/diHMM).

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
