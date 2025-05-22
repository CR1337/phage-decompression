import os
from tqdm import tqdm
from subprocess import Popen
from pathlib import Path
from itertools import product
from multiprocessing import Pool, cpu_count


DIRECTORY = Path("GenomesDB")

fasta_files = [f for f in DIRECTORY.iterdir() if f.is_file() and f.name.endswith(".fasta")]
fasta_files = sorted(fasta_files, key=lambda f: f.name)
print(f"{len(fasta_files)=}")
fasta_files = fasta_files[:64]

if os.path.exists("list.txt"):
    os.remove("list.txt")

def do_decompression(fasta_file):
    name = os.path.basename(fasta_file.name).split(".")[0]
    Popen(
        [
            "python3", "phagetools.py", 
            "decompress", 
            "-i", str(fasta_file),
            "-g", os.path.join(DIRECTORY, name, f"{name}.gff"),
            "-l", "list.txt",
            "-o", os.path.join(DIRECTORY, name)
        ]
    ).wait()

    Popen(
        [
            "python3", "phagetools.py",
            "glm2",
            "-i", str(fasta_file),
            "-g", os.path.join(DIRECTORY, name, f"{name}.gff"),
            "-r",
            "-o", os.path.join(DIRECTORY, name, f"{name}_glm2.txt")
        ]
    ).wait()

    Popen(
        [
            "python3", "phagetools.py",
            "glm2",
            "-i", os.path.join(DIRECTORY, name, f"{name}_decompressed.fasta"),
            "-g", os.path.join(DIRECTORY, name, f"{name}_decompressed.gff"),
            "-o", os.path.join(DIRECTORY, name, f"{name}_glm2_decompressed.txt")
        ]
    ).wait()


with Pool(cpu_count() // 2) as pool:
    list(tqdm(
        pool.imap(do_decompression, fasta_files), 
        total=len(fasta_files), 
        desc="Decompression"
    ))

print("Summary...")
Popen(
    [
        "python3", "phagetools.py",
        "summary",
        "-l", "list.txt",
        "-o", "summary"
    ]
).wait()


def do_similarity(args):
    filename_1, filename_2 = args
    Popen(
        [
            "python3", "phagetools.py",
            "similarity",
            "-i", filename_1, filename_2,
            "-o", "similarities.csv"
        ]
    ).wait()   


with Pool(cpu_count() // 2) as pool:
    list(tqdm(
        pool.imap(do_similarity, product(fasta_files, repeat=2)), 
        total=len(fasta_files) ** 2, 
        desc="Similarity"
    ))


print("Plot similarities...")
Popen(
    [
        "python3", "phagetools.py",
        "similarity",
        "-i", "similarities.csv",
        "-p",
        "-o", "similarities.png"
    ]
).wait()

Popen(
    [
        "python3", "phagetools.py",
        "similarity",
        "-i", "similarities.csv",
        "-p", "-c",
        "-o", "similarities_cluster.png"
    ]
).wait()

print("Done")
