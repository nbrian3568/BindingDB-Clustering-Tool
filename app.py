from flask import Flask, request, render_template, jsonify, send_file, abort
import os
import uuid
from faster_clustering_method import process_bindingdb_tsv

app = Flask(__name__)

UPLOAD_FOLDER = "uploads"
PROCESSED_FOLDER = "processed"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(PROCESSED_FOLDER, exist_ok=True)

app.config["UPLOAD_FOLDER"] = UPLOAD_FOLDER
app.config["PROCESSED_FOLDER"] = PROCESSED_FOLDER
app.config["MAX_CONTENT_LENGTH"] = 1000 * 1024 * 1024  # Max 500MB

@app.route("/", methods=["GET"])
def index():
    return render_template("index.html")

@app.route("/upload", methods=["POST"])
def upload_file():
    if "file" not in request.files:
        return jsonify({"error": "No file part"}), 400
    file = request.files["file"]
    if file.filename == "":
        return jsonify({"error": "No selected file"}), 400
    if not file.filename.endswith(".tsv"):
        return jsonify({"error": "Please upload a valid .tsv file"}), 400

    filename = f"{uuid.uuid4().hex}.tsv"
    filepath = os.path.join(app.config["UPLOAD_FOLDER"], filename)
    file.save(filepath)

    try:
        uniprot_ids = process_bindingdb_tsv(filepath, app.config["PROCESSED_FOLDER"])
    except ValueError as ve:
        return jsonify({"error": str(ve)}), 400
    except Exception as e:
        return jsonify({"error": f"Unexpected server error: {str(e)}"}), 500

    if not uniprot_ids:
        return jsonify({"error": "No valid targets found in the uploaded file."}), 400

    results = []

    for uniprot_id in uniprot_ids:
        target_dir = os.path.join(app.config["PROCESSED_FOLDER"], uniprot_id)
        if not os.path.exists(target_dir):
            continue

        clusters = []
        cluster_files = [f for f in os.listdir(target_dir) if f.endswith(".json") and f.startswith("cluster")]
        cluster_files.sort()

        for cluster_json_file in cluster_files:
            try:
                cluster_id = int(cluster_json_file[len("cluster"):-len(".json")])
            except ValueError:
                continue

            hist_file = f"cluster{cluster_id}_hist.png"
            mol_file = f"cluster{cluster_id}_mol.png"

            clusters.append({
                "cluster_id": cluster_id,
                "json_file": f"/processed/{uniprot_id}/{cluster_json_file}",
                "histogram": f"/processed/{uniprot_id}/{hist_file}" if os.path.exists(os.path.join(target_dir, hist_file)) else None,
                "molecule_image": f"/processed/{uniprot_id}/{mol_file}" if os.path.exists(os.path.join(target_dir, mol_file)) else None
            })

        results.append({
            "uniprot_id": uniprot_id,
            "clusters": clusters
        })

    return jsonify({"results": results})

@app.route("/processed/<uniprot_id>/<filename>")
def serve_processed_file(uniprot_id, filename):
    filepath = os.path.join(app.config["PROCESSED_FOLDER"], uniprot_id, filename)
    if os.path.isfile(filepath):
        return send_file(filepath)
    else:
        abort(404)

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port, debug=False)
