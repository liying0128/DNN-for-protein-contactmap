import autogluon as ag
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.PDB import Structure
from PIL import Image

def download_sequences(protein_id):
    url = "https://www.ncbi.nlm.nih.gov/protein/?term=" + protein_id
    response = requests.get(url)
    response.raise_for_status()
    soup = BeautifulSoup(response.content, "html.parser")
    sequences = []
    for record in soup.find_all("div", class_="sequence-record"):
        sequences.append(record.find("pre").text)
    return sequences


def align_sequences(sequences):
    alignment = MSA(sequences)
    return alignment


def calculate_contact_map(alignment):
    model = ag.models.DeepModel(
        num_classes=2,
        dropout_rate=0.5,
        learning_rate=0.001,
        epochs=100,
    )
    features = []
    labels = []
    for i in range(len(alignment)):
        for j in range(i + 1, len(alignment)):
            features.append([alignment[i][j], alignment[j][i]])
            labels.append(alignment[i][j] == alignment[j][i])
    model.fit(features, labels)
    predictions = model.predict(features)
    return predictions


def save_image(predictions, filename):
    plt.imshow(predictions, cmap="hot")
    plt.savefig(filename)


def handle_exceptional_protein_sequences(sequences):
    sequences = [sequence for sequence in sequences if len(sequence) >= 100]
    for sequence in sequences:
        for amino_acid in sequence:
            if amino_acid not in standard_amino_acids:
                sequences.remove(sequence)
    return sequences


def calculate_weights(alignment):
    weights = []
    for i in range(len(alignment)):
        for j in range(i + 1, len(alignment)):
            weights.append(
                np.exp(
                    -1
                    * (
                        sum(
                            [
                                alignment[i][k] == alignment[j][k]
                                for k in range(len(alignment[0]))
                            ]
                        )
                        / len(alignment[0])
                    )
                )
            )
    return weights


def main():
    protein_id = "P00519"
    num_sequences = 1000

    sequences = download_sequences(protein_id)

    sequences = handle_exceptional_protein_sequences(sequences)

    alignment = align_sequences(sequences)

    weights = calculate_weights(alignment)

    predictions = calculate_contact_map(alignment, weights)

    save_image(predictions, "contact_map.png")


if __name__ == "__main__":
    main()
