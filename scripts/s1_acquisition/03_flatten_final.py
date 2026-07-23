import os
import shutil

# Directory where the flattened files are located after the first flattening step
target_dir = "data/raw_counts"

print("--- Starting final flattening ---")

# os.walk(..., topdown=False) allows us to traverse the directory tree from the bottom up, which is useful for removing directories after processing their files.
for root, dirs, files in os.walk(target_dir, topdown=False):
    for file in files:
        # If the file is already in the root of raw_counts (target_dir), do nothing
        if root == target_dir:
            continue

        # Path to the file
        src_file = os.path.join(root, file)
        
        # Create a new name: Folder_Filename, to avoid duplicates
        parent_folder = os.path.basename(root)
        new_filename = f"{parent_folder}_{file}"
        dst_file = os.path.join(target_dir, new_filename)

        # Move the file
        shutil.move(src_file, dst_file)
        print(f"Flattened: {file} -> {new_filename}")

    # After moving files, try to remove the folder
    # os.rmdir will only remove it if it's empty
    if root != target_dir:
        try:
            os.rmdir(root)
            print(f"Removed folder: {root}")
        except OSError:
            print(f"Folder {root} is not empty, skipping.")

print("--- Done! All files are now at the same level in data/raw_counts ---")