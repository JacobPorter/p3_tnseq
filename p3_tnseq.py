#!/usr/bin/env 
import os, sys
import argparse
import requests
import json
import subprocess

def get_genome(parameters):
    target_file = os.path.join(parameters["output_path"],parameters["gid"]+".fna")
    if not os.path.exists(target_file):
        genome_url= "data_url/genome_sequence/?eq(genome_id,gid)&limit(25000)".replace("data_url",parameters["data_url"]).replace("gid",parameters["gid"])
        print genome_url
        headers = {"accept":"application/sralign+dna+fasta"}
        #print "switch THE HEADER BACK!"
        #headers = {'Content-Type': 'application/x-www-form-urlencoded; charset=utf-8'}
        req = requests.Request('GET', genome_url, headers=headers)
        prepared = req.prepare()
        #pretty_print_POST(prepared)
        s = requests.Session()
        response=s.send(prepared)
        handle = open(target_file, 'w')
        if not response.ok:
            sys.stderr.write("API not responding. Please try again later.\n")
            sys.exit(2)
        else:
            for block in response.iter_content(1024):
                handle.write(block)
    return target_file

def get_annotation(parameters):
    target_file =os.path.join(parameters["output_path"],parameters["gid"]+".gff")
    if not os.path.exists(target_file):
        annotation_url= "data_url/genome_feature/?and(eq(genome_id,gid),eq(annotation,PATRIC),or(eq(feature_type,CDS),eq(feature_type,tRNA),eq(feature_type,rRNA)))&limit(25000)".replace("data_url",parameters["data_url"]).replace("gid",parameters["gid"])
        print annotation_url
        headers = {"accept":"application/cufflinks+gff"}
        #print "switch THE HEADER BACK!"
        #headers = {'Content-Type': 'application/x-www-form-urlencoded; charset=utf-8'}
        req = requests.Request('GET', annotation_url, headers=headers)
        prepared = req.prepare()
        #pretty_print_POST(prepared)
        s = requests.Session()
        response=s.send(prepared)
        handle = open(target_file, 'w')
        if not response.ok:
            sys.stderr.write("API not responding. Please try again later.\n")
            sys.exit(2)
        else:
            for block in response.iter_content(1024):
                handle.write(block)
    return target_file

def get_files(job_data, server_data):
    genome_dirs=[job_data["output_path"]]
    job_data["data_url"]=server_data["data_url"]
    job_data["gid"]=job_data["reference_genome_id"]
    for g in [job_data["reference_genome_id"]]:
        get_annotation(job_data)
        get_genome(job_data)
    return genome_dirs

#contrasts are always defined as either [control,treatment] or [control]
#run transit per contrast
#when allow UI definition of conditions and contrasts will need to clean up names for file friendliness
def run_transit(genome_list, library_dict, parameters):
    contrasts=parameters["contrasts"]
    output_path=parameters["output_path"]
    for genome in genome_list:
        cmd=["transit", parameters["recipe"]]
        for contrast in contrasts:
            output_file=os.path.join(output_path, "_".join([parameters["recipe"]]+contrast)+"_transit.txt")
            cur_cmd=list(cmd) #make a copy
            control_files=[]
            exp_files=[]
            condition=contrast[0]
            for r in library_dict[condition]["replicates"]:
                control_files.append(r[genome["genome"]]["wig"])
            if len(contrast) == 2:
                condition= contrast[1]
                for r in library_dict[condition]["replicates"]:
                    exp_files.append(r[genome["genome"]]["wig"])
            if len(control_files) > 0:
                cur_cmd.append(",".join(control_files))
            else:
                sys.stderr.write("Missing control files for "+parameters["recipe"])
                sys.exit(2)
            if parameters["recipe"]=="resampling":
                if len(exp_files) > 0:
                    cur_cmd.append(",".join(exp_files))
                else:
                    sys.stderr.write("Missing exp files for "+parameters["recipe"])
                    sys.exit(2)
            cur_cmd.append(genome["annotation"])
            cur_cmd.append(output_file)
            print " ".join(cur_cmd)
            subprocess.check_call(cur_cmd) #call transit
        


def run_alignment(genome_list, library_dict, parameters): 
    #modifies library_dict sub replicates to include 'bowtie' dict recording output files
    output_path=parameters["output_path"]
    key_handle=open(os.path.join(parameters["output_path"],"output_keys.txt"), 'w')
    for genome in genome_list:
        genome_link=os.path.join(output_path, os.path.basename(genome["genome"]))
        final_cleanup=[]
        if not os.path.exists(genome_link):
            subprocess.check_call(["ln","-s",genome["genome"],genome_link])
        cmd=["tpp", "-bwa", "bwa", "-ref", genome["genome"]]
        #thread_count=multiprocessing.cpu_count()
        #cmd+=["-p",str(thread_count)]
        #if genome["dir"].endswith('/'):
        #    genome["dir"]=genome["dir"][:-1]
        #genome["dir"]=os.path.abspath(genome["dir"])
        #genome["output"]=os.path.join(output_path,os.path.basename(genome["dir"]))
        for library in library_dict:
            rcount=0
            for r in library_dict[library]["replicates"]:
                cur_cleanup=[]
                rcount+=1
                result_name=library+str(rcount)
                target_dir=output_path
            #    target_dir=os.path.join(genome["output"],library,"replicate"+str(rcount))
            #    target_dir=os.path.abspath(target_dir)
            #    subprocess.call(["mkdir","-p",target_dir])
                cur_cmd=list(cmd)
                if "read2" in r:
                    read_link1=os.path.join(output_path, os.path.basename(r["read1"]))
                    read_link2=os.path.join(output_path, os.path.basename(r["read2"]))
                    if not os.path.exists(read_link1):
                        subprocess.check_call(["ln","-s",r["read1"],read_link1])
                    if not os.path.exists(read_link2):
                        subprocess.check_call(["ln","-s",r["read2"],read_link2])
                    cur_cmd+=["-reads1",read_link1,"-reads2",read_link2]
                    name1=os.path.splitext(os.path.basename(r["read1"]))[0]
                    name2=os.path.splitext(os.path.basename(r["read2"]))[0]
                    key_handle.write("\t".join([name1,name2,result_name])+"\n")
                    base_name=os.path.join(target_dir,result_name)
                else:
                    read_link1=os.path.join(output_path, os.path.basename(r["read1"]))
                    if not os.path.exists(read_link1):
                        subprocess.check_call(["ln","-s",r["read1"],read_link1])
                    cur_cmd+=["-reads1",read_link1]
                    name1=os.path.splitext(os.path.basename(r["read1"]))[0]
                    key_handle.write("\t".join([name1,result_name])+"\n")
                    base_name=os.path.join(target_dir,result_name)
                sam_file = base_name+".sam"
                wig_file = base_name+".wig"
                cur_cleanup.append(sam_file)
                bam_file=sam_file[:-4]+".bam"
                r[genome["genome"]]={}
                r[genome["genome"]]["bam"]=bam_file
                r[genome["genome"]]["wig"]=wig_file
                cur_cmd+=["-output",base_name]
                if os.path.exists(bam_file):
                    sys.stderr.write(bam_file+" alignments file already exists. skipping\n")
                else:
                    print " ".join(cur_cmd)
                    subprocess.check_call(cur_cmd) #call bowtie2
                if not os.path.exists(bam_file):
                    subprocess.check_call("samtools view -Su "+sam_file+" | samtools sort -o - - > "+bam_file, shell=True)#convert to bam
                    subprocess.check_call("samtools index "+bam_file, shell=True)
                    #subprocess.check_call('samtools view -S -b %s > %s' % (sam_file, bam_file+".tmp"), shell=True)
                    #subprocess.check_call('samtools sort %s %s' % (bam_file+".tmp", bam_file), shell=True)
                for garbage in cur_cleanup:
                    if os.path.exists(garbage):
                        subprocess.call(["rm", garbage])
        key_handle.close()
        for garbage in final_cleanup:
            subprocess.call(["rm", garbage])


def main(server_setup, job_data):
    required_data=["experimental_conditions","read_files","reference_genome_id", "recipe", "contrasts"]
    for data in required_data:
        if not data in job_data or len(job_data[data]) == 0:
            sys.stderr.write("Missing "+ data +"\n")
            sys.exit(2)
    
    library_dict={}
    library_list=[]
    library_list=job_data["experimental_conditions"]
    output_path=job_data["output_path"]=os.path.abspath(job_data["output_path"])
    #for lib in library_list:
    #    library_dict[lib]={"library":lib}
    #job_data["read_files"]=job_data["read_files"].split()
    #count=0
    #add read/replicate structure to library dict
    #for read in job_data["read_files"]:
    #    replicates=read.split(',')
    #    rep_store=library_dict[library_list[count]]["replicates"]=[]
    #    for rep in replicates:
    #        pair=rep.split('%')
    #        pair_dict={"read1":pair[0]}
    #        if len(pair) == 2:
    #            pair_dict["read2"]=pair[1]
    #        rep_store.append(pair_dict)
    #    count+=1
   
    library_dict = job_data["read_files"] 
    genome_dirs=get_files(job_data, server_setup)
    genome_list=[]
    for g in genome_dirs:
        cur_genome={"genome":[],"annotation":[],"dir":g,"hisat_index":[]}
        for f in os.listdir(g):
            if f.endswith(".fna") or f.endswith(".fa") or f.endswith(".fasta"):
                cur_genome["genome"].append(os.path.abspath(os.path.join(g,f)))
            elif f.endswith(".gff"):
                cur_genome["annotation"].append(os.path.abspath(os.path.join(g,f)))

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
        #if args.index:
        #    if len(cur_genome["hisat_index"]) != 1:
        #        sys.stderr.write("Missing hisat index tar file for "+g+"\n")
        #        sys.exit(2)
        #    else:
        #        cur_genome["hisat_index"]=cur_genome["hisat_index"][0]


        genome_list.append(cur_genome)
    output_path=os.path.abspath(output_path)
    if not os.path.exists(output_path):
        subprocess.call(["mkdir","-p",output_path])
    run_alignment(genome_list, library_dict, job_data)
    run_transit(genome_list, library_dict, job_data)
    #cleanup(genome_list, library_dict, parameters, output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    jobinfo = parser.add_mutually_exclusive_group(required=True)
    jobinfo.add_argument('--jfile', 
            help='json file for job {"reference_genome_id":[x],"experimental_conditions":[x], "transit_params":{key:value}, "output_path":x, "read_files":x')
    jobinfo.add_argument('--jstring', help='json string from user input')
    serverinfo = parser.add_mutually_exclusive_group(required=True)
    serverinfo.add_argument('--sfile', help='server setup JSON file')
    serverinfo.add_argument('--sstring', help='server setup JSON string')
    parser.add_argument('-o', help='output directory. defaults to current directory.', required=False)

    #parser.add_argument('-g', help='genome ids to get *.fna and annotation *.gff', required=True)
    #parser.add_argument('-L', help='csv list of library names for comparison', required=True)
    #parser.add_argument('-p', help='JSON formatted parameter list for TRANSIT keyed to program', required=True)
    #parser.add_argument('-o', help='output directory. defaults to current directory.', required=False)
    #parser.add_argument('readfiles', nargs='+', help="whitespace sep list of read files. shoudld be \
    #        ws separates control (first) from experiment files (second),\
    #        a comma separates replicates, and a percent separates pairs.")
    if len(sys.argv) ==1:
        parser.print_help()
        sys.exit(2)
    args = parser.parse_args()
    try:
        job_data = json.loads(args.jstring) if args.jstring else json.load(open(args.jfile,'r'))
    except:
        sys.stderr.write("Failed to parse user provided form data \n")
        raise
    #parse setup data
    try:
        server_setup= json.loads(args.sstring) if args.sstring else json.load(open(args.sfile,'r'))
    except:
        sys.stderr.write("Failed to parse server data\n")
        raise
    job_data["output_path"]=args.o 
    main(server_setup, job_data)
