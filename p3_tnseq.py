
#!/usr/bin/env 
import os, sys
import argparse
import requests


def get_genome(genome_id):


def get_annotation(parameters):
    annotation_url= "data_url/genome_feature/?and(eq(genome_id,gid),eq(annotation,PATRIC),or(eq(feature_type,CDS),eq(feature_type,tRNA),eq(feature_type,rRNA)))&limit(25000)".replace("data_url",parameters["data_url"]).replace("gid",parameters["gid"])
    headers = {"accept":"application/cufflinks+gff"}
    #print "switch THE HEADER BACK!"
    #headers = {'Content-Type': 'application/x-www-form-urlencoded; charset=utf-8'}
    req = requests.Request('GET', annotation_url, headers=headers)
    prepared = req.prepare()
    #pretty_print_POST(prepared)
    s = requests.Session()
    response=s.send(prepared)
    if not response.ok:
        sys.stderr.write("API not responding. Please try again later.\n")
        sys.exit(2)


def main(genome_list, library_dict, parameters_file, output_dir, gene_matrix=False):

if __name__ == "__main__":
    #modelling input parameters after rockhopper
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', help='csv list of directories each containing a genome file *.fna and annotation *.gff', required=True)
    parser.add_argument('-L', help='csv list of library names for comparison', required=False)
    parser.add_argument('-p', help='JSON formatted parameter list for tuxedo suite keyed to program', required=False)
    parser.add_argument('-o', help='output directory. defaults to current directory.', required=False)
    #parser.add_argument('-x', action="store_true", help='run the gene matrix conversion and create a patric expression object', required=False)
    parser.add_argument('readfiles', nargs='+', help="whitespace sep list of read files. shoudld be \
            in corresponding order as library list. ws separates libraries,\
            a comma separates replicates, and a percent separates pairs.")
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    library_dict={}
    library_list=[]
    if args.L:
        library_list=args.L.strip().split(',')
    #create library dict
    if not len(library_list): library_list.append("results")
    gene_matrix=True
    if not args.o:
        output_dir="./"
    else:
        output_dir=args.o
    for lib in library_list:
        library_dict[lib]={"library":lib}
    count=0
    #add read/replicate structure to library dict
    for read in args.readfiles:
        replicates=read.split(',')
        rep_store=library_dict[library_list[count]]["replicates"]=[]
        for rep in replicates:
            pair=rep.split('%')
            pair_dict={"read1":pair[0]}
            if len(pair) == 2:
                pair_dict["read2"]=pair[1]
            rep_store.append(pair_dict)
        count+=1
    genome_dirs=args.g.strip().split(',')
    genome_list=[]
    for g in genome_dirs:
        cur_genome={"genome":[],"annotation":[],"dir":g,"hisat_index":[]}
        for f in os.listdir(g):
            if f.endswith(".fna") or f.endswith(".fa") or f.endswith(".fasta"):
                cur_genome["genome"].append(os.path.abspath(os.path.join(g,f)))
            elif f.endswith(".gff"):
                cur_genome["annotation"].append(os.path.abspath(os.path.join(g,f)))
            elif f.endswith(".ht2.tar"):
                cur_genome["hisat_index"].append(os.path.abspath(os.path.join(g,f)))

        if len(cur_genome["genome"]) != 1:
            sys.stderr.write("Too many or too few fasta files present in "+g+"\n")
            sys.exit(2)
        else:
            cur_genome["genome"]=cur_genome["genome"][0]
        if len(cur_genome["annotation"]) != 1:
            sys.stderr.write("Too many or too few gff files present in "+g+"\n")
            sys.exit(2)
        else:
            cur_genome["annotation"]=cur_genome["annotation"][0]
        if args.index:
            if len(cur_genome["hisat_index"]) != 1:
                sys.stderr.write("Missing hisat index tar file for "+g+"\n")
                sys.exit(2)
            else:
                cur_genome["hisat_index"]=cur_genome["hisat_index"][0]


        genome_list.append(cur_genome)
    
    main(genome_list,library_dict,args.p,output_dir,gene_matrix)
