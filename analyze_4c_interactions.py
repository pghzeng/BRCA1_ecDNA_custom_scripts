"""Count 4C-like targeted interactions from pairtools .pairs.gz files.

Parameters:
- --vp: BED file defining viewpoint regions. Each non-comment line must contain at least
    four columns: chrom, start, end, name.
- --target: BED file defining target regions. Each non-comment line must contain at least
    four columns: chrom, start, end, name.
- --config: whitespace-delimited sample manifest. Each non-comment line must contain
    two fields: sample_name and path_to_pairs_gz.
- --output: output TSV path. The script writes one row per sample and VP-target combination
    with columns: Sample, Interaction_Name, Interactions, Total_UU, RPM.

Counting logic:
- Only pairtools records with pair_type == UU are counted.
- A read pair contributes to a VP-target interaction if one end falls in the VP interval
    and the other end falls in the target interval, regardless of read order.
- RPM is computed as interaction_count / total_UU * 1,000,000 for each sample.
"""

import gzip
import argparse

def read_bed(file_path):
    regions = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            # Expected BED columns: chrom, start, end, name.
            regions.append({
                'chr': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'name': parts[3]
            })
    return regions

def main():
    parser = argparse.ArgumentParser(description='4C-seq Targeted Interaction Counter')
    parser.add_argument('--vp', required=True, help='Path to the viewpoint BED file')
    parser.add_argument('--target', required=True, help='Path to the target BED file')
    parser.add_argument('--config', required=True, help='Path to the sample manifest file')
    parser.add_argument('--output', required=True, help='Path to the output TSV file')
    args = parser.parse_args()

    # 1. Load region definitions.
    vps = read_bed(args.vp)
    targets = read_bed(args.target)
    
    # Precompute counters for all VP-target combinations.
    combinations = []
    for v in vps:
        for t in targets:
            combinations.append({'vp': v, 'target': t, 'name': f"{v['name']}_{t['name']}"})

    # 2. Read the sample manifest and process samples.
    with open(args.output, 'w') as out_f:
        out_f.write("Sample\tInteraction_Name\tInteractions\tTotal_UU\tRPM\n")
        
        with open(args.config, 'r') as config_f:
            for line in config_f:
                if line.startswith('#') or not line.strip():
                    continue
                sample_name, file_path = line.strip().split()
                
                print(f">>> Processing Sample: {sample_name}")
                
                total_uu = 0
                counts = {c['name']: 0 for c in combinations}
                
                # 3. Read the .pairs.gz file.
                with gzip.open(file_path, 'rt') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        
                        cols = line.split('\t')
                        # pair_type is the 8th column (index 7).
                        if cols[7] != 'UU':
                            continue
                        
                        total_uu += 1
                        
                        # Parse coordinates as (chr1, pos1, chr2, pos2).
                        c1, p1, c2, p2 = cols[1], int(cols[2]), cols[3], int(cols[4])
                        
                        # Check whether this pair matches any VP-target combination.
                        for comb in combinations:
                            v = comb['vp']
                            t = comb['target']
                            
                            # Logic: (R1 in VP and R2 in target) OR (R1 in target and R2 in VP).
                            match = (
                                (c1 == v['chr'] and v['start'] <= p1 <= v['end'] and 
                                 c2 == t['chr'] and t['start'] <= p2 <= t['end'])
                                or
                                (c1 == t['chr'] and t['start'] <= p1 <= t['end'] and 
                                 c2 == v['chr'] and v['start'] <= p2 <= v['end'])
                            )
                            
                            if match:
                                counts[comb['name']] += 1
                
                # 4. Compute and write the results.
                for comb in combinations:
                    c_name = comb['name']
                    inter_count = counts[c_name]
                    rpm = (inter_count / total_uu * 1000000) if total_uu > 0 else 0
                    out_f.write(f"{sample_name}\t{c_name}\t{inter_count}\t{total_uu}\t{rpm:.4f}\n")
                
                print(f"    Done. Total UU: {total_uu}")

if __name__ == "__main__":
    main()