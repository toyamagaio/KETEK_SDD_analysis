import os

base_path = "/mnt/c/Data/"

for folder_name in os.listdir(base_path):
    if '\r' in folder_name:
        clean_name = folder_name.replace('\r', '')
        old_path = os.path.join(base_path, folder_name)
        new_path = os.path.join(base_path, clean_name)
        os.rename(old_path, new_path)
        print(f"Renamed: {folder_name} -> {clean_name}")
