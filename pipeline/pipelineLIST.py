from pathlib import Path
import sys
import csv

parent_dir = Path(__file__).resolve().parents[1]
here = Path(__file__).resolve().parent
if str(parent_dir) not in sys.path:
	sys.path.append(str(parent_dir))

from pipeline.api_intactUNI import intact_uniprot_pipeline

from pipeline.isolate_interactors import isolate_all_uniprot_ids
from pipeline.anotate_interactors import annotate_pipeline, count_interpro


def isolate_ids(list_families, delimiter=';', col_number=0):
	out_path = here / "uniprot_ids_list"
	out_path.mkdir(parents=True, exist_ok=True)
	for file in list_families:
		file_stem = Path(file).stem
		print(file_stem)
		outfile = f'{file_stem}_uniprot.txt'
		uniprot_ids =[]
		with open(here/ file, 'r', encoding='utf-8') as b83:
			reader = csv.reader(b83, delimiter =delimiter)
			for row in reader:
				if row and not row[0].startswith('#'):
					if row[col_number]:
						uniprot_ids.append(row[col_number])

		with open(out_path/outfile, 'w', encoding='utf-8') as out:
			for uniID in uniprot_ids:
				out.write(uniID + '\n')

def partners(group_name, intact_dir, out_path):
    """
	TODO: napisz
    """
    print(f"Starting collection for group: {group_name}")

    intact_dir = Path(intact_dir)
    out_path = Path(out_path)
    out_path.mkdir(parents=True, exist_ok=True)
    out_file = out_path / f"{group_name}_partners.csv"
    
    all_uniprot_ids, partner_pub_dict = isolate_all_uniprot_ids(
        group_name, intact_dir, add_publication_col=True
    )
    
    group_uniprot_ids = set()
    for f in intact_dir.glob(f"{group_name}_uniprot.txt"):
        with open(f, "r") as fin:
            group_uniprot_ids.update(line.strip() for line in fin if line.strip())

    partners = all_uniprot_ids - group_uniprot_ids

    print(f"Writing structured output to {out_file}")
    
    with open(out_file, "w") as f:
        f.write("PartnerUniProtID,Publications\n")
        
        for pid in sorted(partners):
            pub_set = partner_pub_dict.get(pid, set())
            pubs_string = "|".join(sorted(pub_set))
            f.write(f"{pid},{pubs_string}\n")

    print(f"Partner data written successfully.")

def pipeline():
	'''
	Główna funkcja uruchamiająca cały pipeline przetwarzania danych
	Przed uruchomieniem wprowadzić w linii:
	- 80 pliki są poprawne
	- 84 jest prawidłowy separator oraz kolumna z UniProt ids jest dobrze zdefiniowana
	'''
	data_path = parent_dir / "data"
	
	# 1. specify files
	file_names = ["data.tsv"]

	# 2. isolate the names 
	list_families = [data_path / file_name for file_name in file_names]
	isolate_ids(list_families, delimiter='\t', col_number=0)
	
	# 4. Fetch IntAct interaction data per UniProt ID
	print('fetching intact')
	intact_uniprot_pipeline(
		ids_dir=here / "uniprot_ids_list",
		out_dir=here / 'intact_list',
		file_pattern = '_uniprot'
	)
	print('finished fetching')

	# 5. isolate interactors
	for fam in list_families:
		file_stem = Path(fam).stem
		partners(file_stem, intact_dir=here / 'intact_list', out_path=here/'isolated_partners/LIST')

	# 6. uniprot --> interpro domains
	annotate_pipeline(
		input_dir=here/'isolated_partners/LIST'
	)

	# 7. count the domains
	count_interpro(
		input_dir=here/'isolated_partners/LIST', 
		output_dir=here/'isolated_partners/LIST'
	)

if __name__ == '__main__':
	pipeline()