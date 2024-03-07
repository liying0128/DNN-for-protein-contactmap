Introduction
This program is used to extract protein structure information from PDB files and calculate protein contact maps.

Usage
Install the dependent libraries:
pip install -r requirements.txt
Run the program:
python main.py
Parameters
protein_id: Protein ID, defaults to P00519.
num_sequences: Number of sequences to download, defaults to 1000.
Output
The program will generate a file named contact_map.png, which contains an image of the protein contact map.

Notes
This program uses the Biopython library to parse PDB files.
This program uses the MSA library to align protein sequences.
This program uses the AutoGluon library to train the DNN model.
