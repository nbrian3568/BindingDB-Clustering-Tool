import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, Draw
from rdkit.DataStructs import TanimotoSimilarity
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
import json
import os

morgan_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def process_bindingdb_tsv(file_path, output_base_dir):
    df = pd.read_csv(file_path, sep="\t", on_bad_lines="skip")
    if len(df) < 2:
        raise ValueError("File contains fewer than 2 ligands")

    processed_uniprot_ids = []

    # Group by Target Name
    target_groupings = df.groupby("Target Name")
    valid_targets = {id_ for id_, group in target_groupings if len(group) >= 2}
    if not valid_targets:
        raise ValueError("No targets with ≥2 ligands found")

    for target_name, group in target_groupings:
        if len(group) < 2:
            continue

        uniprot_ids = group["UniProt (SwissProt) Primary ID of Target Chain"].dropna().unique()
        if len(uniprot_ids) == 0:
            continue
        uniprot_id = uniprot_ids[0]
        processed_uniprot_ids.append(uniprot_id)

        target_dir = os.path.join(output_base_dir, uniprot_id)
        os.makedirs(target_dir, exist_ok=True)

        ligand_smiles = group["Ligand SMILES"].dropna().tolist()
        affinities = group["Kd (nM)"].dropna().tolist()
        if not ligand_smiles:
            continue

        parsed = [(i, Chem.MolFromSmiles(sm)) for i, sm in enumerate(ligand_smiles)]
        valid = [(i, mol) for i, mol in parsed if mol is not None]

        valid_indices = [i for i, mol in valid]
        molecules = [mol for i, mol in valid]
        smiles_list = [ligand_smiles[i] for i in valid_indices]

        if not molecules or len(molecules) < 2:
            continue

        affinities = [affinities[i] for i in valid_indices if i < len(affinities)]

        fps = [morgan_generator.GetFingerprint(mol) for mol in molecules]
        n = len(molecules)

        distance_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                sim = TanimotoSimilarity(fps[i], fps[j])
                dist = 1 - sim
                distance_matrix[i, j] = distance_matrix[j, i] = dist

        if distance_matrix.size == 0 or distance_matrix.shape[0] < 2:
            raise ValueError("Not enough ligands to perform clustering")

        condensed_dist = squareform(distance_matrix)
        if condensed_dist.size == 0:
            raise ValueError("Condensed distance matrix is empty")

        linkage_matrix = sch.linkage(condensed_dist, method="average")
        clusters = fcluster(linkage_matrix, t=0.4, criterion="distance")

        cluster_affinities = {}
        for idx, cluster_id in enumerate(clusters):
            if idx < len(affinities):
                try:
                    affinity_val = float(affinities[idx])
                    if affinity_val > 0:
                        log_aff = np.log10(affinity_val)
                        cluster_affinities.setdefault(str(cluster_id), []).append(log_aff)
                except (ValueError, TypeError):
                    continue

        cluster_mols = {}
        cluster_smiles = {}
        for idx, cluster_id in enumerate(clusters):
            mol = molecules[idx]
            sm = smiles_list[idx]
            cluster_mols.setdefault(str(cluster_id), []).append(mol)
            cluster_smiles.setdefault(str(cluster_id), []).append(sm)

        cluster_data = {
            "target_name": target_name,
            "uniprot_id": uniprot_id,
            "clusters": []
        }

        for cluster_id, mols in cluster_mols.items():
            if len(mols) < 2:
                continue
            cluster_id_int = int(cluster_id)
            smiles = cluster_smiles[cluster_id]

            img = Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200, 200), legends=smiles)
            mol_image_filename = f"cluster{cluster_id}_mol.png"
            mol_image_path = os.path.join(target_dir, mol_image_filename)
            img.save(mol_image_path)

            aff_list = cluster_affinities.get(cluster_id, [])
            hist_filename = f"cluster{cluster_id}_hist.png"
            hist_path = os.path.join(target_dir, hist_filename)
            plt.figure(figsize=(10, 10))
            if aff_list:
                plt.hist(aff_list, bins=20, color="blue", edgecolor="black")
                plt.title(f"Cluster {cluster_id} Affinity Histogram")
                plt.xlabel("Affinity (log Kd in nM)")
                plt.ylabel("Frequency")
            else:
                plt.text(1.0, 1.0, "No affinity data", ha='center', va='center')
                plt.axis('off')
            plt.savefig(hist_path)
            plt.close()

            cluster_json = {
                "cluster_id": cluster_id_int,
                "affinities": aff_list,
                "count": len(mols),
                "image": mol_image_filename,
                "histogram": hist_filename,
                "smiles": smiles
            }

            cluster_json_path = os.path.join(target_dir, f"cluster{cluster_id}.json")
            with open(cluster_json_path, "w") as cj:
                json.dump(cluster_json, cj, indent=2)

            cluster_data["clusters"].append(cluster_json)

        cluster_data_path = os.path.join(target_dir, "clusterData.json")
        with open(cluster_data_path, "w") as f:
            json.dump(cluster_data, f, indent=2)

        print(f"Saved cluster data for {uniprot_id} to {cluster_data_path}")

    return processed_uniprot_ids
