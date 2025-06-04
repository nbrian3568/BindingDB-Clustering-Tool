document.addEventListener('DOMContentLoaded', () => {
  const form = document.getElementById('uploadForm');
  const statusDiv = document.getElementById('status');
  const resultsDiv = document.getElementById('results');
  const dropdownContainer = document.getElementById('uniprotDropdownContainer');

  let uniprotData = [];  // Store results for all UniProt IDs

  form.addEventListener('submit', async (event) => {
    event.preventDefault();
    console.log("Form submitted");

    const fileInput = document.getElementById('tsvFile');
    if (fileInput.files.length === 0) {
      alert('Please select a file!');
      return;
    }

    const file = fileInput.files[0];
    const formData = new FormData();
    formData.append('file', file);

    statusDiv.textContent = 'Processing file, please wait...';
    resultsDiv.innerHTML = '';
    dropdownContainer.innerHTML = '';  // Clear dropdown

    try {
      const response = await fetch('/upload', {
        method: 'POST',
        body: formData
      });

      if (!response.ok) {
        const errText = await response.text();
        throw new Error('Server error: ' + errText);
      }

      const data = await response.json();

      statusDiv.textContent = 'Processing complete!';

      if (data.error) {
        resultsDiv.textContent = 'Error: ' + data.error;
        return;
      }

      // Save results data for later use
      uniprotData = data.results;

      // Create dropdown
      if (uniprotData.length === 0) {
        statusDiv.textContent = 'No UniProt data found in file.';
        return;
      }

      const label = document.createElement('label');
      label.textContent = 'Select UniProt ID: ';
      label.setAttribute('for', 'uniprotDropdown');

      const select = document.createElement('select');
      select.id = 'uniprotDropdown';

      // Add default option
      const defaultOption = document.createElement('option');
      defaultOption.value = '';
      defaultOption.textContent = '-- Choose UniProt ID --';
      select.appendChild(defaultOption);

      // Populate dropdown options
      uniprotData.forEach(result => {
        const option = document.createElement('option');
        option.value = result.uniprot_id;
        option.textContent = result.uniprot_id;
        select.appendChild(option);
      });

      dropdownContainer.appendChild(label);
      dropdownContainer.appendChild(select);

      // Show results for selected UniProt ID
      select.addEventListener('change', () => {
        const selectedId = select.value;
        resultsDiv.innerHTML = '';  // Clear previous results

        if (!selectedId) {
          // No selection, clear results
          return;
        }

        // Find the data for selected UniProt ID
        const selectedResult = uniprotData.find(r => r.uniprot_id === selectedId);
        if (!selectedResult) {
          resultsDiv.textContent = 'No data found for selected UniProt ID.';
          return;
        }

        // Render clusters for selected UniProt ID
        const uniprotContainer = document.createElement('div');
        uniprotContainer.className = 'uniprot-block';

        const uniprotTitle = document.createElement('h2');
        uniprotTitle.textContent = `UniProt ID: ${selectedResult.uniprot_id}`;
        uniprotContainer.appendChild(uniprotTitle);

        for (const cluster of selectedResult.clusters) {
          const clusterContainer = document.createElement('div');
          clusterContainer.className = 'cluster';

          const clusterTitle = document.createElement('h3');
          clusterTitle.textContent = `Cluster ${cluster.cluster_id}`;
          clusterContainer.appendChild(clusterTitle);

          if (cluster.affinities && cluster.affinities.length > 0) {
            const affinityText = document.createElement('p');
            affinityText.innerHTML = `<strong>Affinities:</strong> ${cluster.affinities.map(a => a.toFixed(2)).join(', ')}`;
            clusterContainer.appendChild(affinityText);
          }

          if (cluster.molecule_image) {
            const molImg = document.createElement('img');
            molImg.src = cluster.molecule_image;
            molImg.alt = 'Molecule Image';
            molImg.style.maxWidth = '400px';
            clusterContainer.appendChild(molImg);
          }

          if (cluster.histogram) {
            const histImg = document.createElement('img');
            histImg.src = cluster.histogram;
            histImg.alt = 'Cluster Histogram';
            histImg.style.maxWidth = '400px';
            clusterContainer.appendChild(histImg);
          }

          if (cluster.json_file) {
            const jsonLink = document.createElement('a');
            jsonLink.href = cluster.json_file;
            jsonLink.textContent = 'Download Cluster JSON';
            jsonLink.target = '_blank';
            clusterContainer.appendChild(jsonLink);
          }

          uniprotContainer.appendChild(clusterContainer);
        }

        resultsDiv.appendChild(uniprotContainer);
      });

    } catch (error) {
      statusDiv.textContent = '';
      resultsDiv.textContent = 'Upload failed: ' + error.message;
    }
  });
});
