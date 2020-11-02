# diHMM v1.0
Here presents an improved implementation of diHMM, contributed by Stephanos Tsoucas and Yan Kai.

diHMM is a novel chromatin segmentation for annotating chromatin states at two length scales. It was originally developed and described in the publication *Multi-scale chromatin state annotation using a hierarchical hidden Markov model* ([Marco et.al 2017](https://www.nature.com/articles/ncomms15011)). Detailed information of the project can be found at [this page](https://github.com/gcyuan/diHMM).

We applied the improved diHMM algorithm here to generate the multi-scale chromatin state annotations for the 127 human reference epigenomes in the [Roadmap and ENCODE consortia](http://www.roadmapepigenomics.org). Detailed information of the information on the 127 epigenomes can be found at [this site](https://egg2.wustl.edu/roadmap/web_portal/meta.html).

## Accessing the multi-scale chromatin state maps in the 127 epigenomes
We generated the chromatin state maps at the nucleosome (200bp resolution) and domain (4kb resolution) level. These maps can be freely downloaded from [here](https://www.dropbox.com/sh/85nxvu1hiwhwm9r/AAB0pQFvwD1KRqpwOOHf6A_Xa?dl=0).

| Epigenome   ID (EID) | Nucleosome | Domain   |        GROUP     |        Standardized Epigenome name                               |        ANATOMY     |
|----------------------|------------|----------|------------------|------------------------------------------------------------------|--------------------|
| E017                 | [download](https://www.dropbox.com/s/uwnio5ibi0639dn/E017_nucleosome.bed.gz?dl=0)   | [download](https://www.dropbox.com/s/pzmr9qiij7rkytl/E017_domain.bed.gz?dl=0) | IMR90            | IMR90   fetal lung fibroblasts Cell Line                         | LUNG               |
| E002                 | download   | download | ESC              | ES-WA7   Cells                                                   | ESC                |
| E008                 | download   | download | ESC              | H9   Cells                                                       | ESC                |
| E001                 | download   | download | ESC              | ES-I3   Cells                                                    | ESC                |
| E015                 | download   | download | ESC              | HUES6   Cells                                                    | ESC                |
| E014                 | download   | download | ESC              | HUES48   Cells                                                   | ESC                |
| E016                 | download   | download | ESC              | HUES64   Cells                                                   | ESC                |
| E003                 | download   | download | ESC              | H1   Cells                                                       | ESC                |
| E024                 | download   | download | ESC              | ES-UCSF4   Cells                                                 | ESC                |
| E020                 | download   | download | iPSC             | iPS-20b   Cells                                                  | IPSC               |
| E019                 | download   | download | iPSC             | iPS-18   Cells                                                   | IPSC               |
| E018                 | download   | download | iPSC             | iPS-15b   Cells                                                  | IPSC               |
| E021                 | download   | download | iPSC             | iPS   DF 6.9 Cells                                               | IPSC               |
| E022                 | download   | download | iPSC             | iPS   DF 19.11 Cells                                             | IPSC               |
| E007                 | download   | download | ES-deriv         | H1   Derived Neuronal Progenitor Cultured Cells                  | ESC_DERIVED        |
| E009                 | download   | download | ES-deriv         | H9   Derived Neuronal Progenitor Cultured Cells                  | ESC_DERIVED        |
| E010                 | download   | download | ES-deriv         | H9   Derived Neuron Cultured Cells                               | ESC_DERIVED        |
| E013                 | download   | download | ES-deriv         | hESC   Derived CD56+ Mesoderm Cultured Cells                     | ESC_DERIVED        |
| E012                 | download   | download | ES-deriv         | hESC   Derived CD56+ Ectoderm Cultured Cells                     | ESC_DERIVED        |
| E011                 | download   | download | ES-deriv         | hESC   Derived CD184+ Endoderm Cultured Cells                    | ESC_DERIVED        |
| E004                 | download   | download | ES-deriv         | H1   BMP4 Derived Mesendoderm Cultured Cells                     | ESC_DERIVED        |
| E005                 | download   | download | ES-deriv         | H1   BMP4 Derived Trophoblast Cultured Cells                     | ESC_DERIVED        |
| E006                 | download   | download | ES-deriv         | H1   Derived Mesenchymal Stem Cells                              | ESC_DERIVED        |
| E062                 | download   | download | Blood   & T-cell | Primary   mononuclear cells from peripheral blood                | BLOOD              |
| E034                 | download   | download | Blood   & T-cell | Primary   T cells from peripheral blood                          | BLOOD              |
| E045                 | download   | download | Blood   & T-cell | Primary   T cells effector/memory enriched from peripheral blood | BLOOD              |
| E033                 | download   | download | Blood   & T-cell | Primary   T cells from cord blood                                | BLOOD              |
| E044                 | download   | download | Blood   & T-cell | Primary   T regulatory cells from peripheral blood               | BLOOD              |
| E043                 | download   | download | Blood   & T-cell | Primary   T helper cells from peripheral blood                   | BLOOD              |
| E039                 | download   | download | Blood   & T-cell | Primary   T helper naive cells from peripheral blood             | BLOOD              |
| E041                 | download   | download | Blood   & T-cell | Primary   T helper cells PMA-I stimulated                        | BLOOD              |
| E042                 | download   | download | Blood   & T-cell | Primary   T helper 17 cells PMA-I stimulated                     | BLOOD              |
| E040                 | download   | download | Blood   & T-cell | Primary   T helper memory cells from peripheral blood 1          | BLOOD              |
| E037                 | download   | download | Blood   & T-cell | Primary   T helper memory cells from peripheral blood 2          | BLOOD              |
| E048                 | download   | download | Blood   & T-cell | Primary   T CD8+ memory cells from peripheral blood              | BLOOD              |
| E038                 | download   | download | Blood   & T-cell | Primary   T helper naive cells from peripheral blood             | BLOOD              |
| E047                 | download   | download | Blood   & T-cell | Primary   T CD8+ naive cells from peripheral blood               | BLOOD              |
| E029                 | download   | download | HSC   & B-cell   | Primary   monocytes from peripheral blood                        | BLOOD              |
| E031                 | download   | download | HSC   & B-cell   | Primary   B cells from cord blood                                | BLOOD              |
| E035                 | download   | download | HSC   & B-cell   | Primary   hematopoietic stem cells                               | BLOOD              |
| E051                 | download   | download | HSC   & B-cell   | Primary   hematopoietic stem cells G-CSF-mobilized Male          | BLOOD              |
| E050                 | download   | download | HSC   & B-cell   | Primary   hematopoietic stem cells G-CSF-mobilized Female        | BLOOD              |
| E036                 | download   | download | HSC   & B-cell   | Primary   hematopoietic stem cells short term culture            | BLOOD              |
| E032                 | download   | download | HSC   & B-cell   | Primary   B cells from peripheral blood                          | BLOOD              |
| E046                 | download   | download | HSC   & B-cell   | Primary   Natural Killer cells from peripheral blood             | BLOOD              |
| E030                 | download   | download | HSC   & B-cell   | Primary   neutrophils from peripheral blood                      | BLOOD              |
| E026                 | download   | download | Mesench          | Bone   Marrow Derived Cultured Mesenchymal Stem Cells            | STROMAL_CONNECTIVE |
| E049                 | download   | download | Mesench          | Mesenchymal   Stem Cell Derived Chondrocyte Cultured Cells       | STROMAL_CONNECTIVE |
| E025                 | download   | download | Mesench          | Adipose   Derived Mesenchymal Stem Cell Cultured Cells           | FAT                |
| E023                 | download   | download | Mesench          | Mesenchymal   Stem Cell Derived Adipocyte Cultured Cells         | FAT                |
| E052                 | download   | download | Myosat           | Muscle   Satellite Cultured Cells                                | MUSCLE             |
| E055                 | download   | download | Epithelial       | Foreskin   Fibroblast Primary Cells skin01                       | SKIN               |
| E056                 | download   | download | Epithelial       | Foreskin   Fibroblast Primary Cells skin02                       | SKIN               |
| E059                 | download   | download | Epithelial       | Foreskin   Melanocyte Primary Cells skin01                       | SKIN               |
| E061                 | download   | download | Epithelial       | Foreskin   Melanocyte Primary Cells skin03                       | SKIN               |
| E057                 | download   | download | Epithelial       | Foreskin   Keratinocyte Primary Cells skin02                     | SKIN               |
| E058                 | download   | download | Epithelial       | Foreskin   Keratinocyte Primary Cells skin03                     | SKIN               |
| E028                 | download   | download | Epithelial       | Breast   variant Human Mammary Epithelial Cells (vHMEC)          | BREAST             |
| E027                 | download   | download | Epithelial       | Breast   Myoepithelial Primary Cells                             | BREAST             |
| E054                 | download   | download | Neurosph         | Ganglion   Eminence derived primary cultured neurospheres        | BRAIN              |
| E053                 | download   | download | Neurosph         | Cortex   derived primary cultured neurospheres                   | BRAIN              |
| E112                 | download   | download | Thymus           | Thymus                                                           | THYMUS             |
| E093                 | download   | download | Thymus           | Fetal   Thymus                                                   | THYMUS             |
| E071                 | download   | download | Brain            | Brain   Hippocampus Middle                                       | BRAIN              |
| E074                 | download   | download | Brain            | Brain   Substantia Nigra                                         | BRAIN              |
| E068                 | download   | download | Brain            | Brain   Anterior Caudate                                         | BRAIN              |
| E069                 | download   | download | Brain            | Brain   Cingulate Gyrus                                          | BRAIN              |
| E072                 | download   | download | Brain            | Brain   Inferior Temporal Lobe                                   | BRAIN              |
| E067                 | download   | download | Brain            | Brain   Angular Gyrus                                            | BRAIN              |
| E073                 | download   | download | Brain            | Brain_Dorsolateral_Prefrontal_Cortex                             | BRAIN              |
| E070                 | download   | download | Brain            | Brain   Germinal Matrix                                          | BRAIN              |
| E082                 | download   | download | Brain            | Fetal   Brain Female                                             | BRAIN              |
| E081                 | download   | download | Brain            | Fetal   Brain Male                                               | BRAIN              |
| E063                 | download   | download | Adipose          | Adipose   Nuclei                                                 | FAT                |
| E100                 | download   | download | Muscle           | Psoas   Muscle                                                   | MUSCLE             |
| E108                 | download   | download | Muscle           | Skeletal   Muscle Female                                         | MUSCLE             |
| E107                 | download   | download | Muscle           | Skeletal   Muscle Male                                           | MUSCLE             |
| E089                 | download   | download | Muscle           | Fetal   Muscle Trunk                                             | MUSCLE             |
| E090                 | download   | download | Muscle           | Fetal   Muscle Leg                                               | MUSCLE_LEG         |
| E083                 | download   | download | Heart            | Fetal   Heart                                                    | HEART              |
| E104                 | download   | download | Heart            | Right   Atrium                                                   | HEART              |
| E095                 | download   | download | Heart            | Left   Ventricle                                                 | HEART              |
| E105                 | download   | download | Heart            | Right   Ventricle                                                | HEART              |
| E065                 | download   | download | Heart            | Aorta                                                            | VASCULAR           |
| E078                 | download   | download | Sm.   Muscle     | Duodenum   Smooth Muscle                                         | GI_DUODENUM        |
| E076                 | download   | download | Sm.   Muscle     | Colon   Smooth Muscle                                            | GI_COLON           |
| E103                 | download   | download | Sm.   Muscle     | Rectal   Smooth Muscle                                           | GI_RECTUM          |
| E111                 | download   | download | Sm.   Muscle     | Stomach   Smooth Muscle                                          | GI_STOMACH         |
| E092                 | download   | download | Digestive        | Fetal   Stomach                                                  | GI_STOMACH         |
| E085                 | download   | download | Digestive        | Fetal   Intestine Small                                          | GI_INTESTINE       |
| E084                 | download   | download | Digestive        | Fetal   Intestine Large                                          | GI_INTESTINE       |
| E109                 | download   | download | Digestive        | Small   Intestine                                                | GI_INTESTINE       |
| E106                 | download   | download | Digestive        | Sigmoid   Colon                                                  | GI_COLON           |
| E075                 | download   | download | Digestive        | Colonic   Mucosa                                                 | GI_COLON           |
| E101                 | download   | download | Digestive        | Rectal   Mucosa Donor 29                                         | GI_RECTUM          |
| E102                 | download   | download | Digestive        | Rectal   Mucosa Donor 31                                         | GI_RECTUM          |
| E110                 | download   | download | Digestive        | Stomach   Mucosa                                                 | GI_STOMACH         |
| E077                 | download   | download | Digestive        | Duodenum   Mucosa                                                | GI_DUODENUM        |
| E079                 | download   | download | Digestive        | Esophagus                                                        | GI_ESOPHAGUS       |
| E094                 | download   | download | Digestive        | Gastric                                                          | GI_STOMACH         |
| E099                 | download   | download | Other            | Placenta   Amnion                                                | PLACENTA           |
| E086                 | download   | download | Other            | Fetal   Kidney                                                   | KIDNEY             |
| E088                 | download   | download | Other            | Fetal   Lung                                                     | LUNG               |
| E097                 | download   | download | Other            | Ovary                                                            | OVARY              |
| E087                 | download   | download | Other            | Pancreatic   Islets                                              | PANCREAS           |
| E080                 | download   | download | Other            | Fetal   Adrenal Gland                                            | ADRENAL            |
| E091                 | download   | download | Other            | Placenta                                                         | PLACENTA           |
| E066                 | download   | download | Other            | Liver                                                            | LIVER              |
| E098                 | download   | download | Other            | Pancreas                                                         | PANCREAS           |
| E096                 | download   | download | Other            | Lung                                                             | LUNG               |
| E113                 | download   | download | Other            | Spleen                                                           | SPLEEN             |
| E114                 | download   | download | ENCODE2012       | A549   EtOH 0.02pct Lung Carcinoma Cell Line                     | LUNG               |
| E115                 | download   | download | ENCODE2012       | Dnd41   TCell Leukemia Cell Line                                 | BLOOD              |
| E116                 | download   | download | ENCODE2012       | GM12878   Lymphoblastoid Cells                                   | BLOOD              |
| E117                 | download   | download | ENCODE2012       | HeLa-S3   Cervical Carcinoma Cell Line                           | CERVIX             |
| E118                 | download   | download | ENCODE2012       | HepG2   Hepatocellular Carcinoma Cell Line                       | LIVER              |
| E119                 | download   | download | ENCODE2012       | HMEC   Mammary Epithelial Primary Cells                          | BREAST             |
| E120                 | download   | download | ENCODE2012       | HSMM   Skeletal Muscle Myoblasts Cells                           | MUSCLE             |
| E121                 | download   | download | ENCODE2012       | HSMM   cell derived Skeletal Muscle Myotubes Cells               | MUSCLE             |
| E122                 | download   | download | ENCODE2012       | HUVEC   Umbilical Vein Endothelial Primary Cells                 | VASCULAR           |
| E123                 | download   | download | ENCODE2012       | K562   Leukemia Cells                                            | BLOOD              |
| E124                 | download   | download | ENCODE2012       | Monocytes-CD14+   RO01746 Primary Cells                          | BLOOD              |
| E125                 | download   | download | ENCODE2012       | NH-A   Astrocytes Primary Cells                                  | BRAIN              |
| E126                 | download   | download | ENCODE2012       | NHDF-Ad   Adult Dermal Fibroblast Primary Cells                  | SKIN               |
| E127                 | download   | download | ENCODE2012       | NHEK-Epidermal   Keratinocyte Primary Cells                      | SKIN               |
| E128                 | download   | download | ENCODE2012       | NHLF   Lung Fibroblast Primary Cells                             | LUNG               |
| E129                 | download   | download | ENCODE2012       | Osteoblast   Primary Cells                                       | BONE               |










Also, they can be accessed in the following way in terminal:
```
# specify the ID of the epigenome that you're interested in
EID=E001
# specify which level of state map you'd like to download. Could be "nucleosome" or "domain"
Type=domain
wget https://genomebrowser-uploads.hms.harvard.edu/data/yk233/diHMM/annotated/$EID_$Type.bed.gz
```
The information on epigenome ID can be found [here](https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15).

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
