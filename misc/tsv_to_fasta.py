input_path = input("Входной TSV: ")
output_path = input("Куда сохранить FASTA: ")

with open(input_path, "r", encoding="utf-8") as f:
    lines = f.readlines()

header = lines[0].strip().split("\t")
seq_idx = header.index("sequence")
id_idx = header.index("sequence_id") if "sequence_id" in header else None

with open(output_path, "w", encoding="utf-8") as out:
    for line in lines[1:]:
        cols = line.strip().split("\t")
        seq = cols[seq_idx]
        seq_id = cols[id_idx] if id_idx is not None else f"seq_{lines.index(line)}"
        out.write(f">{seq_id}\n{seq}\n")

print("Готово!")
