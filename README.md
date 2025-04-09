# Causal Analysis Experiments

This repository contains the experimental implementation described in the [paper](https://pubs.acs.org/doi/10.1021/acs.molpharmaceut.5c00262), focused on causal analysis using data extracted from this [study](https://pubmed.ncbi.nlm.nih.gov/38496642/).

## Overview

The code performs causal inference analysis based on experimental data, generating results and visualizations consistent with the findings presented in the referenced paper.

## Requirements

Before running the code, make sure all required Python packages are installed. You can do this by running:

```bash
pip install -r requirements.txt
```

## Usage

### Run the Main Analysis

To execute the core analysis:

```bash
python main.py
```

### Generate Causal Heatmaps

To produce heatmaps showing causal relationships across different thresholds (as shown in the paper):

```bash
python plot.py
```

## Notes

- Ensure your Python environment matches the versions specified in `requirements.txt`.
- All figures and outputs will be saved to the appropriate directory specified in the code.
