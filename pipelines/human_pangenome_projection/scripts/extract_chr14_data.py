
import gzip
from pathlib import Path

def extract_fasta(input_path, output_path, target_header_start=">14"):
    print(f"Extracting {target_header_start} from {input_path} to {output_path}...")
    
    opener = gzip.open if str(input_path).endswith('.gz') else open
    mode = 'rt' if str(input_path).endswith('.gz') else 'r'
    
    found = False
    with opener(input_path, mode) as f_in, open(output_path, 'w') as f_out:
        recording = False
        for line in f_in:
            if line.startswith('>'):
                # Check if this is the target chromosome
                # Exact match or starts with target + space to avoid matching 140, 14_etc
                if line.strip() == target_header_start or line.startswith(target_header_start + ' ') or line.startswith(target_header_start + '\t'):
                    recording = True
                    found = True
                    print(f"  Found header: {line.strip()}")
                    f_out.write(line)
                else:
                    recording = False
            elif recording:
                f_out.write(line)
    
    if not found:
        print(f"  WARNING: {target_header_start} not found in {input_path}")
    else:
        print("  Done.")

def extract_gff(input_path, output_path, target_chr="14"):
    print(f"Extracting chr {target_chr} from {input_path} to {output_path}...")
    
    opener = gzip.open if str(input_path).endswith('.gz') else open
    mode = 'rt' if str(input_path).endswith('.gz') else 'r'
    
    count = 0
    with opener(input_path, mode) as f_in, open(output_path, 'w') as f_out:
        # Write headers
        f_out.write("##gff-version 3\n")
        
        for line in f_in:
            if line.startswith('#'):
                # Keep comments/headers? Maybe just top headers.
                # The existing gff3 likely has ##sequence-region which we might want to keep or filter
                # For simplicity, let's keep non-sequence-region headers or just basic ones
                if line.startswith('##sequence-region'):
                    if f" {target_chr} " in line:
                        f_out.write(line)
                else:
                     f_out.write(line)
                continue
            
            parts = line.split('\t')
            if len(parts) > 0 and parts[0] == target_chr:
                f_out.write(line)
                count += 1
                
    print(f"  Extracted {count} features.")

def main():
    root = Path(__file__).resolve().parent.parent
    data_dir = root / "cluster_test_data"
    out_dir = root / "alignment_viewer/backend/data"
    
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. CHM13 Fasta
    extract_fasta(
        data_dir / "CHM13.unmasked.fa",
        out_dir / "CHM13.chr14.fa",
        target_header_start=">14"
    )
    
    # 2. GRCh38 Fasta
    extract_fasta(
        data_dir / "Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        out_dir / "GRCh38.chr14.fa",
        target_header_start=">14"
    )
    
    # 3. CHM13 GFF3
    extract_gff(
        data_dir / "CHM13_XY.gff3",
        out_dir / "CHM13.chr14.gff3",
        target_chr="14"
    )
    
    # 4. GRCh38 GFF3
    extract_gff(
        data_dir / "Homo_sapiens.GRCh38.115.chr.gff3",
        out_dir / "GRCh38.chr14.gff3",
        target_chr="14"
    )

if __name__ == "__main__":
    main()
