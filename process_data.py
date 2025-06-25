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
    try:
        print("Step 1: Reading TSV")
        df = pd.read_csv(file_path, sep="\t", on_bad_lines="skip", low_memory=False)
        print("Step 2: Finished reading TSV")
    except Exception as e:
        print(f"[ERROR] Reading {file_path}: {e}")
        return []

    processed_uniprot_ids = []

    for target_name, group in df.groupby("Target Name"):
        uniprot_ids = group["UniProt (SwissProt) Primary ID of Target Chain"].dropna().unique()
        if len(uniprot_ids) != 1 or "/" in target_name:
            continue
        uniprot_id = uniprot_ids[0]
        target_dir = os.path.join(output_base_dir, uniprot_id)
        os.makedirs(target_dir, exist_ok=True)

        smiles_list = group["Ligand SMILES"].dropna().tolist()
        if not smiles_list:
            continue

        affinity_col = next((col for col in ["Kd (nM)", "Ki (nM)", "IC50 (nM)", "EC50 (nM)"]
                             if col in group.columns and group[col].notna().any()), None)
        if affinity_col is None:
            continue

        parsed = []
        for idx, sm in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(sm)
            if mol is None:
                continue
            raw_aff = group.iloc[idx][affinity_col]
            try:
                if isinstance(raw_aff, str):
                    val = float(raw_aff.strip(">~< "))
                else:
                    val = float(raw_aff)
                if val > 0:
                    parsed.append((idx, mol, sm, val))
            except:
                continue

        if len(parsed) < 2:
            continue

        mols = [p[1] for p in parsed]
        fps = [morgan_generator.GetFingerprint(m) for m in mols]
        n = len(mols)

        dmat = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                sim = TanimotoSimilarity(fps[i], fps[j])
                d = 1 - sim
                dmat[i, j] = dmat[j, i] = d

        try:
            linkage = sch.linkage(squareform(dmat), method="average")
            clusters = fcluster(linkage, t=0.4, criterion="distance")
        except:
            clusters = np.ones(n, dtype=int)

        # ✅ Build cluster_dict first
        cluster_dict = {}
        for i, cid in enumerate(clusters):
            cluster_dict.setdefault(cid, []).append(parsed[i])

        # ✅ Now calculate total ligand count from valid clusters
        total_ligands = sum(len(items) for items in cluster_dict.values() if len(items) >= 2)

        cluster_data = {
            "target_name": target_name,
            "uniprot_id": uniprot_id,
            "affinity_type": affinity_col,
            "ligand_count": total_ligands,  # ✅ added here
            "clusters": []
        }

        for cid, items in cluster_dict.items():
            if len(items) < 2:
                continue
            smiles = [i[2] for i in items]
            affinities = [np.log10(i[3]) for i in items]
            mols_in_cluster = [i[1] for i in items]

            # Histogram
            hist_file = f"cluster{cid}_hist.png"
            plt.figure()
            plt.hist(affinities, bins=np.linspace(min(affinities), max(affinities), 20), color="skyblue", edgecolor="black")
            plt.title(f"Cluster {cid} Affinity Distribution ({uniprot_id})")
            plt.xlabel(f"Affinity (log10 {affinity_col})")
            plt.ylabel("Frequency")
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(os.path.join(target_dir, hist_file))
            plt.close()

            # Representative: highest affinity = lowest raw value
            best_item = min(items, key=lambda x: x[3])
            rep_file = f"cluster{cid}_rep.png"
            try:
                rep_img = Draw.MolToImage(best_item[1], size=(200, 200), legend=best_item[2])
                rep_img.save(os.path.join(target_dir, rep_file))
            except:
                rep_file = None

            # Molecule grid
            grid_file = f"cluster{cid}_mol.png"
            try:
                img = Draw.MolsToGridImage(mols_in_cluster, molsPerRow=4, subImgSize=(200, 200), legends=smiles)
                img.save(os.path.join(target_dir, grid_file))
            except:
                grid_file = None

            cluster_json = {
                "cluster_id": int(cid),
                "affinities": [float(a) for a in affinities],
                "count": int(len(items)),
                "image": grid_file,
                "histogram": hist_file,
                "smiles": smiles,
                "representative": rep_file,
                "representative_affinity": float(np.log10(best_item[3]))
            }

            cluster_data["clusters"].append(cluster_json)

            with open(os.path.join(target_dir, f"cluster{cid}.json"), "w") as f:
                json.dump(cluster_json, f, indent=2)

        with open(os.path.join(target_dir, "clusterData.json"), "w") as f:
            json.dump(cluster_data, f, indent=2)

        processed_uniprot_ids.append(uniprot_id)

    return processed_uniprot_ids
