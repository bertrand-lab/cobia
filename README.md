## *cobia*: Prediction and consequences of cofragmentation in metaproteomics

*cobia* is a computational tool to predict the number of cofragmenting ions in a mass spectrometry experiment, focused on metaproteomics. Specifically, we tried to examine the influence of *co*fragmentation *bia*s on metaproteomics.

### Requirements

The installation instructions assume a linux environment, although they can be adapted to work on a Windows or MacOS system as needed. The following python modules are required:

* `pandas` version > 0.2.0
* `pyteomics`
* Retention time prediction:
  - `BioLCCC` *or* RTModel/RTPredict from OpenMS
* Python 2.7 (for compatability with BiolCCC)

### Installation

To install *cobia*, download the [source code](https://github.com/bertrand-lab/cobia/archive/master.zip), and in this directory run the
`setup.py` file. Running the following command will build and install *cobia*:

`python setup.py install`

I recommend installing on a [`conda`](https://conda.io/docs/) environment:

`python setup.py install --user`

Here is a conda environment command that would setup your environment:

```python
conda create --name pyteo_27 python=2.7 pandas
conda activate pyteo_27
pip install pyteomics.biolccc
pip install pyteomics
```

### Peptide Modifications

Peptide modifications are done using [pyteomics](https://pyteomics.readthedocs.io/en/latest/). The following modifications are fixed an applied to all peptides:

1) Peptides digested with trypsin, with no missed cleavages. 
2) Peptides with unknown amino acids ('\*' or 'X') or containing selenocysteine ('U') are removed.
3) Peptides have fixed chemical modifications. The following modifications are currently in place: oxidation of methionine and carbamidomethylation of cysteine.
4) Charge states are assigned deterministically. If a peptide has an internal histidine or arginine/lysine followed by a proline, it is assigned a charge state of 3. Otherwise it is assigned a charge state of 2.

### Retention Time Prediction

Retention time prediction can be done either with [BioLCCC](http://theorchromo.ru/docs/) (part of [pyteomics](https://pyteomics.readthedocs.io/en/latest/)) or using the OpenMS tool RTModel and RTPredict. If you want to use RTModel or RTPredict, [OpenMS](https://www.openms.de/) should be installed.

BioLCCC retention time prediction requires an LC parameter file containing the following characteristics (but see the BioLCCC documentation for details):

* column_length
* column_diameter
* second_solvent_concentration_a
* second_solvent_concentration_b
* gradient_0
* gradient_1
* gradient_2
* flow_rate
* code_format
* linear
* model

```python
cobia peptide_mod_biolccc_rt_prediction -f potential_proteome.fasta -l lc_parameter_file.csv -g custom_gradient_file.csv -n output_name
```

Alternatively, you can train a support vector machine model with RTModel. This requires a mass spectrometry experiment being completed, which will result in a subset of the peptides in the sample to be observed. These retention times and peptides will be used to train an SVM, and RTPredict can be used to then predict the retention times of the potential proteome. Two examples can be found in `testing-scripts/` (above), called 'svm_models_kleiner.sh'; but the documentation for OpenMS RTModel talks in more detail. 

Once you have trained RTModel, you can digest your potential protome with:

```python
cobia database_trpysin -f potential_proteome.fasta -n output_name
```

You can then predict the retention times with RTPredict. Subsequently, the RTPredict output needs to be modified and mass calculated using:

```python
cobia openms_modelled_rt -rtmodel_output.csv -n output_name
```

This output can then be used for cofragmentation prediction as below. 

### Cofragmention Prediction

We use a relatively simple approach, essentially counting the number of isobaric (+/- precursor selection window / 2) and co-eluting peptides from the *potential*  proteome (ie. the fasta file). One key parameter is the ion peak width. In the paper, we've estimated the ion peak width with a simple linear model derived from Hsieh *et al* (2013) mean ion peak width as a function of gradient and column length (`peak width ~ 0.01978 - 0.0005563*(Column Length) + 0.0065488*(Gradient Length).

#### Global Cofragmentation Prediction

In a global approach, the number of cofragmenting ions for all peptides in the potential proteome (via the fasta file) will be determined. This approach is computationally intensive (depending on the size of your fasta file), and so two approaches are used to increase speed: sparse sampling and parallelization. Sparse sampling approximates the number of cofragmenting ions by subsampling every *n*th injection bin. The number of cores can be set as well. Here is an example command:

```python
cobia cofrag_prediction -f lc_predicted_output.csv -l dda_parameters.csv -output_name --global global
```

#### Targeted Cofragmentation Prediction

In a targeted approach, the number of cofragmenting ions for a subset of peptides in the potential proteome will be determined. This approach is much faster than the global approach, so we do not use sparse sampling or parallelization. Peptides must be supplied as a .csv file with the column header 'pep_seq'. For example:

```python
cobia cofrag_prediction -f lc_predicted_output.csv -l dda_parameters.csv -n output_name --global targeted -t target_peps.csv
```
