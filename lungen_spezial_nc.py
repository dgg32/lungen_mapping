### python lungen_spezial_nc.py spezial_fasta spezial_bam nc_bam

import sys
import pysam
import os
from multiprocessing import Pool, Lock
import multiprocessing
import screed

spezial_bin_fasta_folder = sys.argv[1]

spezial_bam = sys.argv[2]


nc_megahit_fasta = sys.argv[3]

nc_bam = sys.argv[4]


fasta_file = "final.contigs.fa"


nc_contig_list = set()

nc_read_list = set()

temp_nc_list_file = nc_bam.replace(".bam", ".list.txt")
if os.path.isfile(temp_nc_list_file):
    for line in open(temp_nc_list_file, 'r'):
        nc_read_list.add(line.strip())

else:
    with screed.open(nc_megahit_fasta) as seqfile:
        for read in seqfile:

            nc_contig_list.add(read.name.split(" ")[0])

    contig_samfile = pysam.AlignmentFile(nc_bam, "rb")

    for c in nc_contig_list:
        nc_read_list.update(set([read.qname for read in contig_samfile.fetch(contig=c)]))
    
    with open(temp_nc_list_file, 'w+') as output:
        output.write("\n".join(nc_read_list))

def worker(contig):

    contig_samfile = pysam.AlignmentFile(spezial_bam, "rb")

    return [read.qname for read in contig_samfile.fetch(contig=contig)]

num_of_cpu = multiprocessing.cpu_count()

extensions = [".fa", ".fasta"]

for (head, dirs, files) in os.walk(spezial_bin_fasta_folder):
    for file in files:
        for extension in extensions:
            if file.endswith(extension):
                current_file_path = os.path.abspath(os.path.dirname(os.path.join(head, file)))
                with_name = current_file_path + "/"+ file

                contig_list = set()

                with screed.open(with_name) as seqfile:
                    for read in seqfile:
                        contig_list.add(read.name.split(" ")[0])


                spezial_samfile = pysam.AlignmentFile(spezial_bam, "rb")

                contigs = [x for x in spezial_samfile.references if x in contig_list]



                with Pool(num_of_cpu) as P:

                    #P.map(worker, contigs)
                    #pass
                    results = P.map(worker, contigs)

                    P.close()
                    P.join()

                    #print (set(results).intersection(nc_read_list))
                    flatten = set([item for sublist in results for item in sublist])
                    print (f"{file}\t{len(flatten)}\t{len(nc_read_list)}\t{len(flatten.intersection(nc_read_list))}")
    