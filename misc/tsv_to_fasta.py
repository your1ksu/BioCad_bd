import csv
import os

input_path = input("Входной TSV: ")
output_path = input("Куда сохранить FASTA (путь к файлу, например data/BCR.fasta): ")

if os.path.isdir(output_path):
    output_path = os.path.join(output_path, "output.fasta")
    print(f"Указана папка, сохраняю как {output_path}")

with open(input_path, "r", encoding="utf-8") as f:
    reader = csv.reader(f, delimiter="\t")
    header = [h.strip('"') for h in next(reader)]
    seq_idx = header.index("sequence")
    id_idx = header.index("sequence_id")

    with open(output_path, "w", encoding="utf-8") as out:
        for row in reader:
            seq_id = row[id_idx].strip('"')
            seq = row[seq_idx].strip('"')
            out.write(f">{seq_id}\n{seq}\n")

print("Готово!")