# Likelihood Recurrent Geometric Network
This is the reference (TensorFlow) implementation of Likelihood-Recurrent Geometric Network. A manuscrip of this work is available on BioRxiv: https://www.biorxiv.org/content/10.1101/2021.09.03.458873v1). All checkpointed models and datasets are available here: https://iowa-my.sharepoint.com/:f:/g/personal/gqi1_uiowa_edu/EkwQUmKpewREj4IvCfXIzfkBy4kLuPtmyIaWAIUH8Tdp5A?e=1MUgwI. Likelihood-RGN builds on Recurrent Geometric Network (RGN) published by Prof. Mohammed AlQuraishi in 2019 [End-to-end differentiable learning of protein structure](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30076-6), by incorporating atomic displacement factors into the deep learning network's loss function.

## Installation and requirements
Extract all files in the [model](https://github.com/mallory-tollefson/rgn) directory in a single location and use `protling.py`, described further below, to train new models and predict structures. Below are the language requirements and package dependencies:

* Python 2.7
* TensorFlow >= 1.4 (tested up to 1.12)
* setproctitle

## Usage
The [`protling.py`](https://github.com/mallory-tollefson/rgn/blob/master/model/protling.py) script facilities training of and prediction using RGN models. Below are typical use cases. The script also accepts a number of command-line options whose functionality can be queried using the `--help` option.

### Train a new model or continue training an existing model
Likelihood-RGN models are described using a configuration file that controls hyperparameters and architectural choices. For a list of available options and their descriptions, see its [documentation](https://github.com/mallory-tollefson/rgn/blob/master/CONFIG.md). Once a configuration file has been created, along with a suitable dataset (download a ready-made [ProteinNetX](TODO) data set or create a new one from scratch using the [`convert_to_tfrecord.py`](https://github.com/mallory-tollefson/rgn/blob/master/data_processing/convert_to_tfrecord.py) script to first create a [ProteinNet](https://github.com/aqlaboratory/proteinnet) data set, followed by updating the [ProteinNet](https://github.com/aqlaboratory/proteinnet) data set to a ProteinNetX data set using the scripts in the [`bfactor_scripts`](https://github.com/mallory-tollefson/rgn/tree/master/bfactor_scripts) directory). The following directory structure must be created:

```
<baseDirectory>/runs/<runName>/<datasetName>/<configurationFile>
<baseDirectory>/data/<datasetName>/[training,validation,testing]
```

Where the first path points to the configuration file and the second path to the directories containing the training, validation, and possibly test sets. Note that `<runName>` and `<datasetName>` are user-defined variables specified in the configuration file that encode the name of the model and dataset, respectively.

Training of a new model can then be invoked by calling:

```
python protling.py <configurationFilePath> -d <baseDirectory>
```

Download a pre-trained model for an example of a correctly defined directory structure. Note that ProteinNet training sets come in multiple "thinnings" and only one should be used at a time by placing it in the main training directory.

To resume training an existing model, run the command above for a previously trained model with saved checkpoints.

### Predict sequences in ProteinNet TFRecords format using a trained model
To predict the structures of proteins already in ProteinNet `TFRecord` format using an existing model with a saved checkpoint, call:

```
python protling.py <configFilePath> -d <baseDirectory> -p -g0
```

This predicts the structures of the dataset specified in the configuration file. By default only the validation set is predicted, but this can be changed using the `-e` option, e.g. `-e weighted_testing` to predict the test set. The `-g0` option sets the GPU to be used to the one with index 0. If a different GPU is available change the setting appropriately.

### Predict structure of a single new sequence using a trained model
If all you have is a single sequence for which you wish to make a prediction, there are multiple steps that must be performed. First, a PSSM needs to be created by running JackHMMer (or a similar tool) against a sequence database, the resulting PSSM must be combined with the sequence in a ProteinNet record, and the file must be converted to the `TFRecord` format. Predictions can then be made as previously described.

Below is an example of how to do this using the supplied scripts (in [data_processing](https://github.com/mallory-tollefson/rgn/blob/master/data_processing)) and one of the pre-trained models, assumed to be unzipped in `<baseDirectory>`. HMMER must also be installed. The raw sequence databases (`<fastaDatabase>`) used in building PSSMs can be obtained from [here](https://github.com/aqlaboratory/proteinnet/blob/master/docs/raw_data.md). The script below assumes that `<sequenceFile>` only contains a single sequence in the FASTA file format.

```
jackhmmer.sh <sequenceFile> <fastaDatabase>
python convert_to_proteinnet.py <sequenceFile>
python convert_to_tfrecord.py <sequenceFile>.proteinnet <sequenceFile>.tfrecord 42
cp <sequenceFile>.tfrecord <baseDirectory>/data/<datasetName>/testing/
python protling.py <baseDirectory>/runs/<runName>/<datasetName>/<configurationFile> -d <baseDirectory> -p -e weighted_testing -g0
```

The first line searches the supplied database for matches to the supplied sequence and extracts a PSSM out of the results. It will generate multiple new files. These are then used in the second line to construct a text-based ProteinNet file (with 42 entries per evolutionary profile, compatible with the pre-trained RGN models). The third line converts the file to `TFRecords` format, and the fourth line copies the file to the testing directory of a pre-trained model. Finally the fifth line predicts the structure using the pre-trained RGN model. The outputs will be placed in  `<baseDirectory>/runs/<runName>/<datasetName>/<latestIterationNumber>/outputsTesting/` and will be comprised of two files: a `.tertiary` file which contains the atomic coordinates, and `.recurrent_states` file which contains the RGN latent representation of the sequence. The `-g0` option sets the GPU to be used to the one with index 0. If a different GPU is available change the setting appropriately.

## Pre-trained models
Below we make available pre-trained Likelihood-RGN models using the [ProteinNetX](TODO) 12 dataset as checkpointed TF graphs and as raw data.

| [CASP12 X-Ray](https://iowa-my.sharepoint.com/:f:/g/personal/gqi1_uiowa_edu/El02wu5TE-hHtHIG74BEFpEB-eokDAKDEXirQuwY6eFTCQ?e=kRjXXq) | [CASP12 X-Ray Raw Data](https://iowa-my.sharepoint.com/:u:/g/personal/gqi1_uiowa_edu/EcbHUDeS061EskTQZ559ZrgB2B8cETdiYQfqpVnrfxAOlQ?e=cinJGb) | [CASP12 X-Ray+NMR](https://iowa-my.sharepoint.com/:f:/g/personal/gqi1_uiowa_edu/EsgIEbshIn9AkzrPs3uw-mgBhwPZcLcvO3DvTJbubfGqfQ?e=ethGZq) | [CASP12 X-Ray+NMR Raw Data](https://iowa-my.sharepoint.com/:u:/g/personal/gqi1_uiowa_edu/EWydoWuMAl9OvshFLOrhSC0BzBLn9ibVjUDVBg5egOervw?e=QBAQgA) |
| --- | --- | --- | --- |

To train new models from scratch using the same hyperparameter choices as the above models, use the appropriate configuration file from [here](https://github.com/mallory-tollefson/rgn/blob/master/configurations/).

## PyTorch implementation
The reference Likelihood-RGN implementation is currently only available in TensorFlow, however the [OpenProtein](https://github.com/OpenProtein/openprotein) project implements various aspects of Prof. AlQuraishi's RGN model in PyTorch, and [PyTorch-RGN](https://github.com/conradry/pytorch-rgn) is a work-in-progress implementation of the RGN model.

## Reference
Our work is available on BioRxiv: https://www.biorxiv.org/content/10.1101/2021.09.03.458873v1

## Funding
This material is based upon work supported by the NSF (National Science Foundation) Graduate Research Fellowship under Grant No. 000390183 to Mallory Tollefson. Guowei Qi was supported by the Barry Goldwater Foundation and the Iowa Center for Research by Undergraduates. Prof. Michael Schnieders was supported by NIH R01DK110023, NIH R01DC012049, and NSF CHE-1751688.
