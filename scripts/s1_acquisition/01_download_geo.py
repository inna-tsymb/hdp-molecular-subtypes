import GEOparse
import os
import signal

base_dir = "data"
metadata_dir = os.path.join(base_dir, "metadata")
raw_files_dir = os.path.join(base_dir, "raw_counts")

os.makedirs(metadata_dir, exist_ok=True)
os.makedirs(raw_files_dir, exist_ok=True)

# Datasets list (Discovery + Validation)
accession_list = [
    "GSE203507", "GSE148241", "GSE306864", "GSE190971", "GSE186257", "GSE75010", "GSE97320", "GSE154377", "GSE192902" # Discovery and Validation datasets
]

def download_data(gse_id):
    print(f"--- Starting download: {gse_id} ---")
    try:
        
        gse = GEOparse.get_GEO(geo=gse_id, destdir=metadata_dir)
        metadata_path = os.path.join(metadata_dir, f"{gse_id}_metadata.csv")
        gse.phenotype_data.to_csv(metadata_path)
        
        print(f"Downloading additional files for {gse_id}...")
        gse.download_supplementary_files(directory=raw_files_dir)
        
        print(f"--- Completed: {gse_id} ---")
        
    except Exception as e:
        print(f"Error downloading {gse_id}: {e}")
        print("Continuing with the next dataset...")

if __name__ == "__main__":
    for gse_id in accession_list:
        download_data(gse_id)
    print("All data downloaded. Check the 'data' folder.")