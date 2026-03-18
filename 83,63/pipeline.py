from pathlib import Path
import sys
import csv


parent_dir = Path(__file__).resolve().parents[1]
here = Path(__file__).resolve().parent
if str(parent_dir) not in sys.path:
	sys.path.append(str(parent_dir))

from pipeline.api_intactUNI import intact_uniprot_pipeline

from pipeline.isolate_interactors import isolate_all_uniprot_ids, write_partner_ids_for_group
from pipeline.anotate_interactors import annotate_pipeline, count_interpro


def isolate_ids(list_families):
	for file in list_families:
		file_stem = Path(file).stem
		print(file_stem)
		outfile = f'{file_stem}_uniprot.txt'
		uniprot_ids =[]
		with open(here/ file, 'r', encoding='utf-8') as b83:
			reader = csv.reader(b83, delimiter =';')
			header = next(reader)
			for row in reader:
				# print(row)
				if row[0]:
					uniprot_ids.append(row[0])

		with open(here/outfile, 'w', encoding='utf-8') as out:
			for uniID in uniprot_ids:
				out.write(uniID + '\n')


def fetch_intact():
	print('fetching intact')
	intact_uniprot_pipeline(
		ids_dir=here,
		out_dir=here,
		file_pattern = '*_uniprot.txt'
	)
	print('finished fetching')

def partners(group_name, intact_dir, out_path):
	print(f"Starting collection of partner UniProt IDs for group: {group_name}")

	intact_dir = Path(intact_dir)
	out_path = Path(out_path)
	out_path.mkdir(parents=True, exist_ok=True)
	out_file = out_path / f"{group_name}_partners.txt"
	print(f"Collecting all UniProt IDs for group '{group_name}' from files in: {intact_dir}")
	all_uniprot_ids = isolate_all_uniprot_ids(group_name, intact_dir)
	

	group_uniprot_ids = set()
	for f in intact_dir.glob(f"{group_name}_uniprot.txt"):
		with open(f, "r") as fin:
			group_uniprot_ids.update(line.strip() for line in fin if line.strip())
	partners = all_uniprot_ids - group_uniprot_ids

	print(f"Total partner UniProt IDs (excluding already present): {len(partners)}")
	with open(out_file, "w") as f:
		for pid in sorted(partners):
			f.write(pid + "\n")
	print(f"Partner UniProt IDs for group '{group_name}' written to {out_file}")

def pipeline():
	list_families = ['83.csv','63.csv']
	# isolate_ids(list_families)
	# fetch_intact()
	# 6. isolate interactors
	# for fam in list_families:
	# 	file_stem = Path(fam).stem
	# 	partners(file_stem, intact_dir=here, out_path=here/'isolated_partners')

	#7. uniprot --> interpro domains
	# annotate_pipeline(
	# 	uniprot_partner_dir=here/'isolated_partners',
	# 	output_dir=here/'isolated_partners'
	# )

	count_interpro(
		input_dir=here/'isolated_partners', 
		output_dir=here/'isolated_partners'
	)

if __name__ == '__main__':
	pipeline()