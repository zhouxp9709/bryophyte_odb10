from Bio import SeqIO
import subprocess
import os
import argparse
import tempfile
import concurrent.futures
from functools import lru_cache
from collections import defaultdict
import sys

@lru_cache(maxsize=10000)
def calculate_similarity_cached(seq1_str, seq2_str, species1, species2):
    if species1 == species2:
        return None

    temp_dir = '/dev/shm' if os.path.exists('/dev/shm') else None
    
    try:
        with tempfile.NamedTemporaryFile(dir=temp_dir, mode="w", delete=False, suffix=".fasta") as f1, \
             tempfile.NamedTemporaryFile(dir=temp_dir, mode="w", delete=False, suffix=".fasta") as f2:
            
            f1.write(f">seq1\n{seq1_str}\n")
            f2.write(f">seq2\n{seq2_str}\n")
            temp_seq1, temp_seq2 = f1.name, f2.name

        with tempfile.NamedTemporaryFile(dir=temp_dir, mode="w", delete=False, suffix=".txt") as f3:
            needle_output = f3.name

        subprocess.run(
            ["needle", "-asequence", temp_seq1, "-bsequence", temp_seq2,
             "-gapopen", "10", "-gapextend", "0.5", "-outfile", needle_output,
             "-brief", "Y"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )

        similarity = 0.0
        with open(needle_output, "r") as f:
            for line in f:
                if line.startswith("# Identity:"):
                    parts = line.strip().split()
                    identical, total = parts[2].split('/')
                    similarity = int(identical) / int(total)
                    break

    except Exception as e:
        similarity = 0.0
    finally:
        for f in [temp_seq1, temp_seq2, needle_output]:
            try:
                os.remove(f)
            except:
                pass

    return similarity

def process_species_pair(args):
    gene_seq_str, species, gene_id, other_species, other_genes_seqs_str = args
    max_sim = 0.0
    
    if not other_genes_seqs_str:
        return (gene_id, other_species, 0.0)
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        futures = [
            executor.submit(
                calculate_similarity_cached,
                gene_seq_str,
                other_seq_str,
                species,
                other_species
            ) for other_seq_str in other_genes_seqs_str
        ]
        
        for future in concurrent.futures.as_completed(futures):
            try:
                sim = future.result()
                if sim and sim > max_sim:
                    max_sim = sim
            except Exception as e:
                pass
    
    return (gene_id, other_species, max_sim)

def read_fasta(file_path):
    sequences = defaultdict(list)
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            if "|" not in record.id:
                continue
            species, gene_id = record.id.split("|", 1)
            sequences[species].append(record)
    except:
        sys.exit(1)
    
    if not sequences:
        sys.exit(1)
    
    return sequences

def find_best_copies(sequences, num_workers, tmpdir):
    os.environ['TMPDIR'] = tmpdir
    species_list = list(sequences.keys())
    best_copies = {}

    if all(len(genes) == 1 for genes in sequences.values()):
        return {species: genes[0] for species, genes in sequences.items()}

    tasks = []
    for species in species_list:
        for gene in sequences[species]:
            gene_seq_str = str(gene.seq)
            gene_id = gene.id
            for other_species in species_list:
                if other_species == species:
                    continue
                other_genes = [str(g.seq) for g in sequences.get(other_species, [])]
                if other_genes:
                    tasks.append((
                        gene_seq_str,
                        species,
                        gene_id,
                        other_species,
                        other_genes
                    ))

    if tasks:
        similarity_data = defaultdict(dict)
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = {executor.submit(process_species_pair, task): task for task in tasks}
            for future in concurrent.futures.as_completed(futures):
                task = futures[future]
                try:
                    gene_id, other_species, max_sim = future.result()
                    similarity_data[(gene_id, other_species)] = max_sim
                except:
                    pass

        for species in species_list:
            if len(sequences[species]) == 1:
                best_copies[species] = sequences[species][0]
                continue

            max_avg = -1
            best_gene = None
            for gene in sequences[species]:
                total = 0.0
                count = 0
                for other_species in species_list:
                    if other_species == species:
                        continue
                    key = (gene.id, other_species)
                    sim = similarity_data.get(key, 0)
                    if sim > 0:
                        total += sim
                        count += 1
                if count > 0:
                    avg = total / count
                    if avg > max_avg:
                        max_avg = avg
                        best_gene = gene
            
            if best_gene:
                best_copies[species] = best_gene
    
    else:
        for species in species_list:
            best_copies[species] = sequences[species][0]

    return best_copies

def write_output(best_copies, output_file):
    try:
        with open(output_file, "w") as f:
            for species, record in best_copies.items():
                SeqIO.write(record, f, "fasta")
    except:
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("-j", "--workers", type=int, default=os.cpu_count())
    parser.add_argument("-t", "--tmpdir", default="/tmp")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        sys.exit(1)
    
    os.makedirs(args.tmpdir, exist_ok=True)
    
    sequences = read_fasta(args.input)
    best_copies = find_best_copies(sequences, args.workers, args.tmpdir)
    write_output(best_copies, args.output)

if __name__ == "__main__":
    main()
