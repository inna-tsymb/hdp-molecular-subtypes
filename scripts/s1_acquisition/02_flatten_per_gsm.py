import os
import shutil

src_root = "data/raw_files"
dst_root = "data/raw_counts"

os.makedirs(dst_root, exist_ok=True)

print("--- Starting flattening and renaming ---")

for folder in os.listdir(src_root):
    folder_path = os.path.join(src_root, folder)
    
    if os.path.isdir(folder_path):
        # Iterate through files in the folder and move them to the destination with new names
        for filename in os.listdir(folder_path):
            src_file = os.path.join(folder_path, filename)
            
            # Make new filename with folder name as prefix
            new_filename = f"{folder}_{filename}"
            dst_file = os.path.join(dst_root, new_filename)
            
            # Move and rename the file
            shutil.move(src_file, dst_file)
            print(f"Flattened: {filename} -> {new_filename}")
        
        # Delete the now-empty folder
        if not os.listdir(folder_path):
            os.rmdir(folder_path)
            print(f"Removed empty folder: {folder}")

print("--- Done! All files in data/raw_counts and structured ---")