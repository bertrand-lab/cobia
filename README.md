## *cobia*: Prediction and consequences of cofragmentation in metaproteomics

*cobia* is a computational tool to predict the number of cofragmenting ions in a mass spectrometry experiment, focused on metaproteomics. Specifically, we calculate 'cofragmentation scores' which represent identification and quantification bias using mass spectrometry based metaproteomics.

### Requirements

The installation instructions assume a linux environment, although they can be adapted to work on a Windows or MacOS system as needed. The following python modules are required:

* `pandas` version > 0.2.0
* `pyteomics`
* Retention time prediction:
  - `BioLCCC` *or* RTModel/RTPredict from OpenMS
* Python 2.7 (for compatability with BiolCCC)
* `pyopenms` (for subsampling retention times from an idXML file)

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
pip install pyopenms
```

To use *cobia*, please see the vignette above!
