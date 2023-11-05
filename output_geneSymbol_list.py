# Any env with scanpy should be fine
# HB: conda activate py-geneformer
# to use:
# python3.8 output_geneSymbol_list.py --inputDir 

import getopt
import argparse
import os
import scanpy as sc


parser = argparse.ArgumentParser(description='output_geneSymbol_list for each h5ad file under an input dir and write the list into the same input dir.')

parser.add_argument('--inputDir', type=str, default='/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/Benchmark_data_top2000/COVID/', help='Directory of input: default(/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/Benchmark_data_top2000/COVID/)')

args = parser.parse_args()

input_dir = args.inputDir


def output_geneSymbol_list(input_dir):
    for file in os.listdir(input_dir):
        if file.endswith(".h5ad"):
            file_base = file.split(".h5ad")[0]
            ann_ob = sc.read_h5ad(input_dir+file)      
            gene_list = ann_ob.var['gene_name'].tolist()
            file_out_name = file_base+"_genesymbol_list.txt"
            # Open the file in write mode ('w')
            with open(input_dir+file_out_name, 'w') as f:
                # Write each item in the list to the file
                for item in gene_list:
                    f.write("%s\n" % item)

            
#input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/Benchmark_data_top2000/COVID/"
output_geneSymbol_list(input_dir=input_dir)
