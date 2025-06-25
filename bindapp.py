from flask import Flask, render_template, abort, url_for
import os
import json

app = Flask(__name__)

DATA_DIR = "static/clusters"

@app.route("/")
def home():
    targets = []
    for uniprot_id in os.listdir(DATA_DIR):
        cluster_file = os.path.join(DATA_DIR, uniprot_id, "clusterData.json")
        if os.path.exists(cluster_file):
            with open(cluster_file) as f:
                target_data = json.load(f)
                targets.append(target_data)
    return render_template("index.html", targets=targets)

@app.route("/clusters/<uniprot_id>")
def show_clusters(uniprot_id):
    cluster_data_file = os.path.join(DATA_DIR, uniprot_id, "clusterData.json")
    if not os.path.exists(cluster_data_file):
        abort(404)
    with open(cluster_data_file) as f:
        data = json.load(f)
    return render_template("clusters.html", target=data)

@app.route("/clusters/<uniprot_id>/<int:cluster_id>")
def show_cluster_details(uniprot_id, cluster_id):
    cluster_file = os.path.join(DATA_DIR, uniprot_id, f"cluster{cluster_id}.json")
    if not os.path.exists(cluster_file):
        abort(404)
    with open(cluster_file) as f:
        cluster = json.load(f)
    return render_template("cluster_view.html", cluster=cluster, uniprot_id=uniprot_id)

if __name__ == "__main__":
    app.run(debug=True)
