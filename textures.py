import os
from magic_dorye.utils.file_meneger import manage_input_dir

def extract_spin_textures(scf_name_and_path):
    scf_name, path = manage_input_dir(scf_name_and_path)
    prefix = scf_name.split(".")[0]
    output_file = os.path.join(path, f"{prefix}.bands1.spin_{direction}.gnu")

    files = {
            'x': os.path.join(path, f"{prefix}.bands1.dat.1"),
            'y': os.path.join(path, f"{prefix}.bands1.dat.2"),
            'z': os.path.join(path, f"{prefix}.bands1.dat.3"),
        }
    
        for direction, file_path in files.items():
            with open(file_path, "r") as f:
                lines = f.read().strip().splitlines()
    
            bands = []
            current_band = []
    
            for line in lines:
                stripped = line.strip()
                if stripped.startswith("&plot") or stripped == "":
                    continue
    
                parts = stripped.split()
                if len(parts) == 3:
                    # k-point line: start a new band if we have data collected
                    if current_band:
                        bands.append(current_band)
                        current_band = []
                else:
                    # Accumulate band data
                    values = [float(x) for x in parts]
                    current_band.extend(values)
    
            # Don't forget the last band
            if current_band:
                bands.append(current_band)
    
            # Write each band as a column with \n in between
            output_file = os.path.join(path, f"{prefix}.bands1.spin_{direction}.gnu")
            with open(output_file, "w") as out:
                for band in bands:
                    for value in band:
                        out.write(f"{value:.6f}\n")
                    out.write("\n")
    
            print(f"Wrote spin-{direction} bands to: {output_file}")

# Example usage
extract_spin_textures("mos2/qe/soc/spin_textures/goodies/mos2.scf.in")
